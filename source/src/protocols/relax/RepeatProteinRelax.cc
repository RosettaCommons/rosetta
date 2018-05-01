// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RepeatProteinRelax.cc
/// @brief Allows the protein to maintain symmetry during internal relax
/// @author TJ Brunette and Phil Bradley

// Unit headers
#include <protocols/relax/RepeatProteinRelax.hh>
#include <protocols/relax/RepeatProteinRelaxCreator.hh>


// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/minimization_packing/MinMover.hh>
// headers likely needed for movers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <boost/foreach.hpp>


namespace protocols {
namespace relax {
static basic::Tracer TR( "protocols.simple_moves.RepeatProteinRelax" );
using namespace core;
using namespace conformation;
using namespace std;
using utility::vector1;
using ObjexxFCL::format::F;


RepeatProteinRelax::RepeatProteinRelax(){
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// gets the midpoint of the longest loop
///
Size RepeatProteinRelax::get_midpoint_longest_loop(core::pose::Pose const & pose, Size const repeat_length){
	protocols::loops::Loops pose_loops;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	string lastSecStruct = dssp_string.substr(0,1);
	Size startLoop = 0;
	Size endLoop = 0;
	if ( dssp_string.substr(0,1) == "L" ) {
		startLoop = 1;
	}
	for ( core::Size ii =repeat_length-5; ii <= repeat_length*2+5; ++ii ) { //Check the second repeat, pad the range by 5
		if ( dssp_string.substr(ii-1,1) == "L" && lastSecStruct != "L" ) {
			startLoop = ii;
		}
		if ( dssp_string.substr(ii-1,1) != "L" && lastSecStruct == "L" ) {
			endLoop = ii-1;
			if ( (startLoop != 1) && (endLoop!=pose.size()) ) {
				pose_loops.add_loop(startLoop,endLoop);
			}
		}
		lastSecStruct = dssp_string.substr(ii-1,1);
	}
	Size max_loop_size=0;
	Size max_loop_position=0;
	for ( Size ii=1; ii<=pose_loops.size(); ii++ ) {
		Size loop_size =pose_loops[ii].stop()-pose_loops[ii].start();
		if ( loop_size>max_loop_size ) {
			max_loop_size=loop_size;
			max_loop_position=ii;
		}
	}
	Size mid_point_longest_loop = (pose_loops[max_loop_position].stop()-pose_loops[max_loop_position].start())/2+pose_loops[max_loop_position].start();
	if ( mid_point_longest_loop>repeat_length ) {
		mid_point_longest_loop=mid_point_longest_loop-repeat_length;
	}
	return(mid_point_longest_loop);
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// fills in a SymmetryInfo object with the necessary info
///
void
RepeatProteinRelax::setup_repeat_symminfo(
	Size const repeatlen,
	Size const nrepeat,
	Size const base_repeat,
	bool const with_jumping,
	conformation::symmetry::SymmetryInfo & symminfo
)
{

	Size const base_offset( (base_repeat-1)*repeatlen );

	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;
		Size const offset( (i-1)*repeatlen );
		for ( Size j=1; j<= repeatlen; ++j ) {
			symminfo.add_bb_clone ( base_offset+j, offset+j );
			symminfo.add_chi_clone( base_offset+j, offset+j );
		}
	}

	symminfo.num_virtuals( 1 ); // the one at the end...
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( nrepeat*repeatlen+1, 1 );
	if ( !with_jumping ) {
		symminfo.torsion_changes_move_other_monomers( true ); // note -- this signals that we are folding between repeats
	}

	/// repeats prior to base_repeat have score_multiply ==> 0
	for ( Size i=1; i<base_repeat; ++i ) {
		Size const offset( (i-1)*repeatlen );
		for ( Size j=1; j<= repeatlen; ++j ) symminfo.set_score_multiply( offset+j, 0 );
	}

	symminfo.update_score_multiply_factor();
}


///////////////////////////////////////////////////////////////////////////////
/// sets up a repeat pose, starting from a non-symmetric pdb with nres=repeatlen*nrepeat
///
void
RepeatProteinRelax::setup_repeat_pose(
	Pose & pose
)
{
	runtime_assert( !pose::symmetry::is_symmetric( pose ) ); // not to begin with...
	Size repeatlen = pose.size()/numb_repeats_;
	Size nrepeat = numb_repeats_;
	Size base_repeat = 2;

	runtime_assert( nrepeat * repeatlen == pose.size() );

	runtime_assert( base_repeat > 1 ); // why? well, the base repeat should get the right context info from nbring repeats
	// but note that with base_repeat>1 we probably can't use linmem_ig and there are other places in the code that
	// assume that monomer 1 is the independent monomer. These should gradually be fixed. Let me (PB) know if you run into
	// trouble.

	Size const nres_protein( nrepeat * repeatlen );
	remove_upper_terminus_type_from_pose_residue( pose, pose.size() );
	ResidueOP vrtrsd
		( conformation::ResidueFactory::create_residue( *core::pose::get_restype_for_pose( pose, "VRTBB" ) ) );
	pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
	pose.conformation().insert_chain_ending( nres_protein );

	kinematics::FoldTree f( pose.size() );
	f.reorder( pose.size() );
	pose.fold_tree( f );


	conformation::symmetry::SymmetryInfo symminfo;

	setup_repeat_symminfo( repeatlen, nrepeat, base_repeat, false, symminfo ); //false sets the with_jumping variable to false

	/// now make symmetric
	pose::symmetry::make_symmetric_pose( pose, symminfo );


	runtime_assert( pose::symmetry::is_symmetric( pose ) );

	//TJ adding these to Phil's function
	TR.Trace << "symmetrize backbone pose" << endl;
	Size const base_offset( (base_repeat-1)*repeatlen );
	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		pose.set_phi  ( pos, pose.phi  (pos) );
		pose.set_psi  ( pos, pose.psi  (pos) );
		pose.set_omega( pos, pose.omega(pos) );
	}

	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		ResidueOP oldrsd( pose.residue(pos).clone() );
		pose.replace_residue( pos, *oldrsd, false );
	}

}


///////////////////////////////////////////////////////////////////////////////
/// sets up a repeat pose with jumping
///
void
RepeatProteinRelax::setup_repeat_pose_jumping(
	Pose & pose
)
{
	runtime_assert( !pose::symmetry::is_symmetric( pose ) ); // not to begin with...
	Size repeatlen = pose.size()/numb_repeats_;
	Size nrepeat = numb_repeats_;
	Size base_repeat = 2;
	Size anchor1 = 2;
	Size anchor2 = repeatlen-2;
	Size cutpoint = get_midpoint_longest_loop(pose,repeatlen);
	runtime_assert( nrepeat * repeatlen == pose.size() );

	runtime_assert( base_repeat > 1 ); // why? well, the base repeat should get the right context info from nbring repeats
	// but note that with base_repeat>1 we probably can't use linmem_ig and there are other places in the code that
	// assume that monomer 1 is the independent monomer. These should gradually be fixed. Let me (PB) know if you run into
	// trouble.

	Size const nres_protein( nrepeat * repeatlen );
	remove_upper_terminus_type_from_pose_residue( pose, pose.size() );

	for ( Size i=1; i<= nrepeat; ++i ) {
		Size const offset( (i-1)*repeatlen);

		ResidueOP vrtrsd
			( conformation::ResidueFactory::create_residue( *core::pose::get_restype_for_pose( pose, "VRT" ) ) );
		//
		Vector ca1( pose.residue(offset+1).xyz("CA")),
			ca2( pose.residue(offset+2).xyz("CA")) ,
			ca3( pose.residue(offset+3).xyz("CA")),
			orig( ca1 ),
			x( (ca2-ca1).normalized() ),
			z( x.cross( ca3-ca1 ).normalized() ),
			y( z.cross(x));
		vrtrsd->set_xyz("ORIG", orig );
		vrtrsd->set_xyz("X", orig+x );
		vrtrsd->set_xyz("Y", orig+y );
		pose.append_residue_by_jump( *vrtrsd, 1 );
	}
	pose.conformation().insert_chain_ending( nres_protein );

	runtime_assert( pose.size() == nres_protein + nrepeat );

	kinematics::FoldTree f( pose.size() );

	for ( Size i=1; i<= nrepeat; ++i ) {
		Size const my_virtual( nres_protein + i );
		Size const offset( (i-1)*repeatlen);
		f.new_jump( offset+anchor1, my_virtual, offset+cutpoint );
		f.new_jump( offset+anchor2, my_virtual, offset+repeatlen );
		//add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, offset+cutpoint );
		//add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, offset+cutpoint+1 );
		if ( i<nrepeat ) f.new_jump( my_virtual, my_virtual+1, my_virtual );
	}
	f.reorder( nres_protein + base_repeat ); // the virtual for the base repeat
	pose.fold_tree( f );


	conformation::symmetry::SymmetryInfo symminfo;

	setup_repeat_symminfo( repeatlen, nrepeat, base_repeat, true, symminfo ); //false sets the with_jumping variable to false

	/// now make symmetric
	pose::symmetry::make_symmetric_pose( pose, symminfo );


	runtime_assert( pose::symmetry::is_symmetric( pose ) );



	//TJ adding these to Phil's function
	TR.Trace << "symmetrize backbone pose" << endl;
	Size const base_offset( (base_repeat-1)*repeatlen );
	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		pose.set_phi  ( pos, pose.phi  (pos) );
		pose.set_psi  ( pos, pose.psi  (pos) );
		pose.set_omega( pos, pose.omega(pos) );
	}

	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		ResidueOP oldrsd( pose.residue(pos).clone() );
		pose.replace_residue( pos, *oldrsd, false );
	}

}


void
RepeatProteinRelax::relax_pose(
	Pose & pose
)
{
	core::pose::PoseOP pdb_poseOP = pose.clone();
	protocols::relax::FastRelax fastrelax( get_scorefxn(), relax_iterations_ );
	if ( cartesian_ ) {
		fastrelax.cartesian(true);
	}
	kinematics::MoveMapOP movemap =  get_movemap();
	fastrelax.set_movemap( movemap );
	TR.Trace << "start relaxing" << endl;
	fastrelax.apply( pose );
	TR.Trace << "finished relaxing" << endl;
	Real const rmsd =  core::scoring::CA_rmsd( pose, *pdb_poseOP);
	TR.Trace << "rmsd to pdb_pose after relax: " << F(9,3,rmsd)<< endl;
}

void
RepeatProteinRelax::minimize_pose(
	Pose & pose
)
{
	kinematics::MoveMapOP movemap =  get_movemap();
	bool const use_nblist( true ), deriv_check( true ), deriv_check_verbose( false );
	protocols::minimization_packing::MinMoverOP min_mover
		( new protocols::minimization_packing::MinMover(movemap, get_scorefxn(), min_type(), 1e-2,
		use_nblist, deriv_check, deriv_check_verbose ) );
	if ( cartesian_ ) {
		min_mover->cartesian(true);
	}
	min_mover->apply( pose );
}

void RepeatProteinRelax::setup_movemap(
	Pose & pose
){
	core::kinematics::MoveMapOP mm = get_movemap();
	mm->set_chi( true );
	mm->set_bb( true );
	mm->set_jump( true );
	//
	//
	mm->set_bb ( pose.size(), false );
	mm->set_chi( pose.size(), false );
	set_movemap(mm);
}

void RepeatProteinRelax::seal_jumps(
	Pose & pose
){
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).name() == "VRTBB" ) {
			pose.conformation().delete_residue_slow(i);
		}
	}
	for ( Size ii=1; ii<=pose.size()-1; ++ii ) {
		if ( pose.residue( ii ).has_variant_type( chemical::CUTPOINT_LOWER ) && pose.residue( ii+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, ii );
			remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, ii+1 );
		}
	}
	kinematics::FoldTree ft = pose.fold_tree();
	ft.clear();
	ft.add_edge(1,pose.size(),core::kinematics::Edge::PEPTIDE);
	pose.fold_tree(ft);
}

void RepeatProteinRelax::add_residue_labels_back(core::pose::Pose & pose, std::map<core::Size, utility::vector1<std::string> > res_label_map,int symmetry_resid_offset){
	std::map<core::Size, vector1<std::string> >::iterator res_label_map_itr;
	vector1<std::string>::iterator res_label_itr;
	for ( res_label_map_itr=res_label_map.begin(); res_label_map_itr!=res_label_map.end(); ++res_label_map_itr ) {
		Size resid = res_label_map_itr->first+symmetry_resid_offset;
		vector1<std::string> tmp_labels = res_label_map_itr->second;
		for ( res_label_itr=tmp_labels.begin(); res_label_itr!=tmp_labels.end(); ++res_label_itr ) {
			pose.pdb_info()->add_reslabel(resid, *res_label_itr);
		}
	}
}


void RepeatProteinRelax::apply(core::pose::Pose & pose) {
	std::map<core::Size, vector1<std::string> > res_label_map;
	Size tmp_pose_size=pose.size();
	for ( core::Size ii = 1; ii <=tmp_pose_size; ++ii ) {
		vector1<std::string> tmp_labels = pose.pdb_info()->get_reslabels(ii);
		res_label_map.insert(std::pair<Size,vector1<std::string> >(ii,tmp_labels));
	}
	if ( modify_symmetry_and_exit_ && remove_symm_ ) {
		pose::symmetry::make_asymmetric_pose(pose);
		seal_jumps(pose);
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
		add_residue_labels_back(pose,res_label_map,0);
		return;
	}
	if ( loop_cutpoint_mode_ ) {
		setup_repeat_pose_jumping(pose);
	} else {
		setup_repeat_pose(pose);
	}
	setup_movemap(pose);
	if ( modify_symmetry_and_exit_ && !remove_symm_ ) {
		add_residue_labels_back(pose,res_label_map,0);
		return;
	}
	if ( minimize_ ) {
		minimize_pose(pose);
	} else {
		relax_pose(pose);
	}
	pose::symmetry::make_asymmetric_pose(pose);
	seal_jumps(pose);
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
	add_residue_labels_back(pose,res_label_map,0);
}

std::string RepeatProteinRelax::get_name() const { return "RepeatProteinRelax"; }



void RepeatProteinRelax::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	numb_repeats_ = tag->getOption<Size>("numb_repeats");
	loop_cutpoint_mode_ = tag->getOption<bool>("loop_cutpoint_mode",false);

	relax_iterations_ = tag->getOption<Size>("relax_iterations",5);
	minimize_ = tag->getOption<bool>("minimize",false);

	//relax options
	cartesian_  = tag->getOption< bool >( "cartesian", false );
	ramp_down_constraints_ = tag->getOption< bool >( "ramp_down_constraints", false );

	if ( tag->hasOption( "min_type" ) ) {
		min_type(tag->getOption< std::string >( "min_type" ));
	} else {
		min_type("lbfgs_armijo_nonmonotone");
	}

	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, datamap )->clone() );

	if ( loop_cutpoint_mode_ ) {
		//get_scorefxn()->set_weight(core::scoring::cart_bonded, 1.0 );
		//get_scorefxn()->set_weight(core::scoring::cart_bonded_length,0.2);
		//get_scorefxn()->set_weight(core::scoring::cart_bonded_angle,0.5);
		get_scorefxn()->set_weight(core::scoring::linear_chainbreak, 2.0);
	}

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi( true );
	mm->set_bb( true );
	mm->set_jump( true );

	remove_symm_  = tag->getOption< bool >( "remove_symmetry", true );
	modify_symmetry_and_exit_= tag->getOption< bool >( "modify_symmetry_and_exit", false );
}

void RepeatProteinRelax::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute("name", xs_string, "Name of the mover")
		+ XMLSchemaAttribute::required_attribute( "numb_repeats", xsct_non_negative_integer, "number of repeats to output" )
		+ XMLSchemaAttribute::attribute_w_default("minimize", xsct_rosetta_bool, "minimize instead of relax","false")
		+ XMLSchemaAttribute::attribute_w_default("loop_cutpoint_mode", xsct_rosetta_bool, "adds a cutpoint to the longest loop but keeps the junctions rigid","false")
		+ XMLSchemaAttribute::attribute_w_default("relax_iterations", xsct_non_negative_integer, "number of iterations to relax","5")
		+ XMLSchemaAttribute::attribute_w_default("cartesian", xsct_rosetta_bool, "cartesian relax","false")
		+ XMLSchemaAttribute::attribute_w_default("min_type", xs_string, "type of minimization","lbfgs_armijo_nonmonotone")
		+ XMLSchemaAttribute::attribute_w_default("remove_symmetry", xsct_rosetta_bool, "removes symmetry setup after mover","true")
		+ XMLSchemaAttribute::attribute_w_default("modify_symmetry_and_exit", xsct_rosetta_bool, "sets up symmetry and exits,note does not delete symmetry","false")
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function used for relax" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "setups internal repeat symmetry or performs repeat symmetry relax depending on call", attlist );
}


std::string RepeatProteinRelax::mover_name(){
	return "RepeatProteinRelax";
}

std::string RepeatProteinRelaxCreator::keyname() const {
	return RepeatProteinRelax::mover_name();
}

protocols::moves::MoverOP
RepeatProteinRelaxCreator::create_mover() const {
	return protocols::moves::MoverOP( new RepeatProteinRelax );
}

void RepeatProteinRelaxCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RepeatProteinRelax::provide_xml_schema( xsd );
}


} // moves
} // protocols
