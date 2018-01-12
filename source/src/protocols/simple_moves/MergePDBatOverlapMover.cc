// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MergePDBatOverlapMover.cc
/// @brief This class combines two pdbs with a known overlap
/// @author TJ Brunette (tjbrunette@gmail.com)
///

// Unit headers
#include <protocols/simple_moves/MergePDBatOverlapMoverCreator.hh>
#include <protocols/simple_moves/MergePDBatOverlapMover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.simple_moves.MergePDBatOverlapMover" );

#include <utility/tag/Tag.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/types.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Conformation.hh>

#include <core/select/util.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/xyzVector.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using core::pose::Pose;
using utility::vector1;

std::string MergePDBatOverlapMoverCreator::keyname() const{
	return MergePDBatOverlapMover::mover_name();
}

protocols::moves::MoverOP
MergePDBatOverlapMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MergePDBatOverlapMover );
}

void MergePDBatOverlapMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MergePDBatOverlapMover::provide_xml_schema( xsd );
}

MergePDBatOverlapMover::MergePDBatOverlapMover()
: moves::Mover("MergePDBatOverlapMover")
{
}

MergePDBatOverlapMover::MergePDBatOverlapMover(core::scoring::ScoreFunctionOP sfxn)
: moves::Mover("MergePDBatOverlapMover")
{
	sfxn_ = sfxn;
}

moves::MoverOP
MergePDBatOverlapMover::clone() const
{
	return moves::MoverOP( new MergePDBatOverlapMover( *this ) );
}

moves::MoverOP
MergePDBatOverlapMover::fresh_instance() const
{
	return moves::MoverOP( new MergePDBatOverlapMover );
}


/// @brief increases the range to ignore to include the entire secondary structure element
void MergePDBatOverlapMover::increase_range_to_ignore_ss_element(core::pose::Pose const & pose, Size init_start, Size init_end, Size & ss_start, Size & ss_end){
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char ss_type = dssp.get_dssp_secstruct(init_start);
	ss_start = init_start;
	while ( ss_type == dssp.get_dssp_secstruct(ss_start) && ss_start>1 ) {
		ss_start--;
	}
	ss_start++;
	ss_end=init_end;
	ss_type = dssp.get_dssp_secstruct( init_end);
	while ( ss_type == dssp.get_dssp_secstruct(ss_end) && ss_end<pose.total_residue() ) {
		ss_end++;
	}
	ss_end--;
}

Size MergePDBatOverlapMover::closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid){
	Real min_dist = 999;
	Size closest_residue = 999999;
	for ( Size ii=1; ii<=pose.total_residue(); ++ii ) {
		if ( ii<start_overlap_resid || ii>end_overlap_resid ) {
			Real tmp_dist = pose.residue(ii).xyz("CA").distance(pose.residue(resid).xyz("CA"));
			if ( tmp_dist < min_dist ) {
				min_dist = tmp_dist;
				closest_residue = ii;
			}
		}
	}
	return(closest_residue);
}



void MergePDBatOverlapMover::merge_junction_sequence(Pose & pose,std::string pose_junction_seq,std::string attach_pose_junction_seq,Size first_overlap_position){
	using namespace core::scoring;
	using namespace core::optimization;
	using namespace core::pack::task;
	bool postions_modified = false;
	Size last_overlap_position=first_overlap_position+pose_junction_seq.size();
	utility::vector1< bool > overlap_and_neighbors( pose.size(), false);
	for ( int ii=0; ii<(int)attach_pose_junction_seq.length(); ++ii ) {
		if ( pose_junction_seq.at(ii) != attach_pose_junction_seq.at(ii) ) {
			postions_modified=true;
			Size resnum = ii+first_overlap_position;
			overlap_and_neighbors[resnum]=true;
			Size ss_start=0;
			Size ss_end=0;
			increase_range_to_ignore_ss_element(pose, first_overlap_position, last_overlap_position, ss_start, ss_end);
			Size closest_residue = closest_non_overlap_residue(pose,resnum,ss_start,ss_end);
			if ( attachment_termini_ == "n_term" ) {
				if ( closest_residue <first_overlap_position ) {
					assign_seq(pose,attach_pose_junction_seq.at(ii),resnum);
				} else {
					assign_seq(pose,pose_junction_seq.at(ii),resnum);
				}
			}
			if ( attachment_termini_ == "c_term" ) {
				if ( closest_residue <first_overlap_position ) {
					assign_seq(pose,pose_junction_seq.at(ii),resnum);
				} else {
					assign_seq(pose,attach_pose_junction_seq.at(ii),resnum);
				}
			}
		}
	}
	if ( postions_modified ) {
		Size packing_range = 5;
		core::select::fill_neighbor_residues(pose, overlap_and_neighbors, packing_range);
		optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
		kinematics::MoveMap mm_loc;
		mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( false );
		TR << "mutations to relax residue to relax" << std::endl;
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( overlap_and_neighbors[ resnum ] ) {
				TR << resnum <<"+";
				mm_loc.set_chi(resnum, true);
			}
		}
		TR << std::endl;
		minopt.max_iter( 100 );
		AtomTreeMinimizer minimizer;
		minimizer.run( pose, mm_loc, *sfxn_, minopt );
	}
}

bool MergePDBatOverlapMover::merge_poses(Pose & pose,Pose & attach_pose){
	using namespace core::id;
	using namespace core::scoring;
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	std::map<core::Size, core::Size> residue_map;
	core::pose::initialize_atomid_map( atom_map, attach_pose, AtomID::BOGUS_ATOM_ID() );
	std::string pose_junction_seq = "";
	std::string attach_pose_junction_seq = "";
	//get the appropriate residues
	if ( attachment_termini_ == "n_term" ) {
		pose_junction_seq = pose.sequence().substr(0,overlap_length_);
		attach_pose_junction_seq = attach_pose.sequence().substr(attach_pose.total_residue()-overlap_length_,overlap_length_);
	} else {
		pose_junction_seq = pose.sequence().substr(pose.total_residue()-overlap_length_,overlap_length_);
		attach_pose_junction_seq = attach_pose.sequence().substr(0,overlap_length_);
	}
	//superimpose---------------------------------------------------------------
	Size start_overlap_pose = 0;
	Size start_overlap_attach_pose = 0;
	if ( attachment_termini_ == "n_term" ) {
		start_overlap_pose = 1;
		start_overlap_attach_pose = attach_pose.size()-overlap_length_+1;
	}
	if ( attachment_termini_ == "c_term" ) {
		start_overlap_pose = pose.size()-overlap_length_+1;
		start_overlap_attach_pose = 1;
	}
	Size exlude_tail_length = 7;
	for ( Size ii=0; ii<overlap_length_-(exlude_tail_length*2); ++ii ) {
		core::id::AtomID const id1(pose.residue(start_overlap_pose+ii+exlude_tail_length).atom_index("CA"),start_overlap_pose+ii+exlude_tail_length);
		core::id::AtomID const id2(attach_pose.residue(start_overlap_attach_pose+ii+exlude_tail_length).atom_index("CA"), start_overlap_attach_pose+ii+exlude_tail_length);
		atom_map[id2]=id1;
		Size residue_pose =start_overlap_pose+ii+exlude_tail_length;
		Size residue_attach = start_overlap_attach_pose+ii+exlude_tail_length;
		residue_map.insert(std::pair<Size,Size>(residue_attach,residue_pose));
	}
	core::Real rmsd = CA_rmsd(attach_pose,pose,residue_map);
	TR << "Junction rmsd" << rmsd << std::endl;
	superimpose_pose(attach_pose,pose,atom_map);
	//attach_pose.dump_pdb("super_attach_pose.pdb");
	//pose.dump_pdb("super_pose.pdb");
	if ( rmsd>max_overlap_rmsd_ ) {
		TR << "junction failed due to rmsd error of" << rmsd << "> then threshold" <<max_overlap_rmsd_ << std::endl;
		return(false);
	}
	//create the poses that will be the connections--------------------------------------------------
	utility::vector1<Size> pose_positions;
	utility::vector1<Size> attach_pose_positions;
	Size n_term_overlap = overlap_length_/2;
	Size first_overlap_position = 0;
	if ( attachment_termini_== "n_term" ) {
		first_overlap_position = attach_pose.total_residue()-overlap_length_+1;
		for ( Size ii=overlap_length_-n_term_overlap; ii<=pose.size(); ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=1; ii<=attach_pose.size()-n_term_overlap-1; ++ii ) {
			attach_pose_positions.push_back(ii);

		}
	}
	if ( attachment_termini_== "c_term" ) {
		first_overlap_position = pose.total_residue()-overlap_length_+1;
		for ( Size ii=1; ii<pose.size()-n_term_overlap; ++ii ) {
			pose_positions.push_back(ii);

		}
		for ( Size ii=overlap_length_-n_term_overlap; ii<=attach_pose.size(); ++ii ) {
			attach_pose_positions.push_back(ii);
		}
	}
	pose::Pose ref_pose_slice;
	pose::Pose attach_pose_slice;
	pdbslice(ref_pose_slice,pose,pose_positions);
	pdbslice(attach_pose_slice,attach_pose,attach_pose_positions);
	Real tmpScore = sfxn_->score(ref_pose_slice);
	std::cout << "score after ref_pose_slice" << tmpScore << std::endl;
	tmpScore = sfxn_->score(attach_pose_slice);
	std::cout << "score after attach_pose_slice" << tmpScore << std::endl;
	//ref_pose_slice.dump_pdb("super_ref_pose_slice.pdb");
	//attach_pose_slice.dump_pdb("super_attach_slice.pdb");
	if ( attachment_termini_== "n_term" ) {
		remove_upper_terminus_type_from_pose_residue(attach_pose_slice,attach_pose_slice.size());
		remove_lower_terminus_type_from_pose_residue(ref_pose_slice,1);
		for ( Size ii=1; ii<=ref_pose_slice.total_residue(); ++ii ) {
			attach_pose_slice.append_residue_by_bond(ref_pose_slice.residue(ii),false/*ideal bonds*/);
		}
		//append_pose_to_pose(attach_pose_slice,ref_pose_slice,false);
		pose = attach_pose_slice;
	}
	if ( attachment_termini_== "c_term" ) {
		remove_upper_terminus_type_from_pose_residue(ref_pose_slice,ref_pose_slice.size());
		remove_lower_terminus_type_from_pose_residue(attach_pose_slice,1);
		for ( Size ii=1; ii<=attach_pose_slice.total_residue(); ++ii ) {
			ref_pose_slice.append_residue_by_bond(attach_pose_slice.residue(ii),false/*ideal bonds*/);
		}
		//append_pose_to_pose(ref_pose_slice,attach_pose_slice,false);
		pose = ref_pose_slice;
	}
	tmpScore = sfxn_->score(pose);
	std::cout << "score after merging" << tmpScore << std::endl;
	//pose.dump_pdb("super_after_attach.pdb");
	if ( pose_junction_seq != attach_pose_junction_seq ) {
		merge_junction_sequence(pose,pose_junction_seq,attach_pose_junction_seq,first_overlap_position);
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
	tmpScore = sfxn_->score(pose);
	std::cout << "score after merging sequence" << tmpScore << std::endl;
	//minimize core
	if ( minimize_after_overlap_ ) {
		minimize_overlap(pose,first_overlap_position,first_overlap_position+overlap_length_);
	}
	return(true);
}

void MergePDBatOverlapMover::assign_seq(Pose & pose, char residue_type, Size position){
	using namespace core::chemical;
	simple_moves::MutateResidueOP mutation_mover;
	AA my_aa =aa_from_oneletter_code(residue_type);
	mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
		position, //position
		my_aa//residue
		) );
	mutation_mover->apply( pose);
	TR << "Ignore the core.conformation.Residue warnings." << std::endl;
}

void MergePDBatOverlapMover::minimize_overlap(Pose & pose,Size overlap_start,Size overlap_end) {
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::optimization;
	core::Real current_score = sfxn_->score(pose);
	std::cout << "current score begin min overlap " << current_score << std::endl;

	utility::vector1< bool > overlap_and_neighbors( pose.size(), false);
	std::cout << "overlap_start" << overlap_start << "overlap_end" << overlap_end << std::endl;
	for ( Size ii=overlap_start; ii<=overlap_end; ++ii ) {
		overlap_and_neighbors[ii]=true;
	}
	Size midpoint=overlap_start+(overlap_end-overlap_start)/2;  // I don't place contraints on the location of the junction so the joint is hopefully minimized
	Size packing_range = 5;
	core::select::fill_neighbor_residues(pose, overlap_and_neighbors, packing_range);
	TR << "minimized residues";
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		if ( overlap_and_neighbors[ resnum ] ) {
			TR << resnum <<"+";
		}
	}
	TR << std::endl;
	//add constraints
	ConstraintCOPs csts;
	sfxn_->set_weight(core::scoring::coordinate_constraint, 0.7);
	sfxn_->set_weight(core::scoring::cart_bonded, 1.0 );  //to repair any cuts.
	sfxn_->set_weight(core::scoring::cart_bonded_length,0.2);
	sfxn_->set_weight(core::scoring::cart_bonded_angle,0.5);
	core::scoring::func::HarmonicFuncOP coord_cst_func = core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 0.0, 0.2 ) );
	core::id::AtomID const anchor_atom( core::id::AtomID( pose.residue(1).atom_index("CA"), 1) );
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		if ( resnum<midpoint-2 || resnum > midpoint+2 ) {
			Size atomindex =  pose.residue(resnum ).atom_index( "CA" );
			csts.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint(core::id::AtomID( atomindex, resnum),anchor_atom,pose.residue(resnum).xyz( "CA" ),coord_cst_func) ));
		} else {
			std::cout << "unconstrained:" << resnum << std::endl;
		}
	}
	current_score = sfxn_->score(pose);
	std::cout << "current score before constraints " << current_score << std::endl;
	pose.add_constraints( csts );
	current_score = sfxn_->score(pose);
	std::cout << "current score before relax" << current_score << std::endl;
	kinematics::FoldTree fold_tree;
	fold_tree = pose.fold_tree();
	std::cout << "fold_tree" << fold_tree << std::endl;
	//minimize loop carefully
	optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
	kinematics::MoveMap mm_loc;
	mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( false );
	TR << "allow residue to move" << std::endl;
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		if ( overlap_and_neighbors[ resnum ] ) { //should only minimize residues that have been allowed to be minimized.
			TR << resnum << "+";
			mm_loc.set_chi( resnum, true);
			mm_loc.set_bb( resnum, true);
		}
	}
	TR << std::endl;
	minopt.max_iter( 200 );
	core::optimization::CartesianMinimizer minimizer;
	minimizer.run( pose, mm_loc, *sfxn_, minopt );
	pose.remove_constraints();
	current_score = sfxn_->score(pose);
	std::cout << "current score" << current_score << std::endl;
	//pose.dump_pdb("super_after_min.pdb");
}


bool MergePDBatOverlapMover::makeJunctions_apply(core::pose::Pose & pose, core::pose::Pose const & attach_pose, Size overlap_length,core::Real max_overlap_rmsd, std::string attachment_termini){
	attach_pose_= attach_pose.clone();
	overlap_length_ = overlap_length;
	max_overlap_rmsd_ = max_overlap_rmsd;
	attachment_termini_ = attachment_termini;
	minimize_after_overlap_ = true;
	return(apply_helper(pose));
}

bool MergePDBatOverlapMover::apply_helper( Pose & pose ){
	bool success = merge_poses(pose,*attach_pose_);
	if ( !success ) {
		return false;
	}
	return true;
}

void MergePDBatOverlapMover::apply( Pose & pose )
{
	bool success = apply_helper(pose );
	if ( !success ) {
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	}
}


void MergePDBatOverlapMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	attachment_termini_ = tag->getOption< std::string >( "attachment_termini" ,"n_term" );
	overlap_length_ = tag->getOption< core::Size >( "overlap_length", 40);
	max_overlap_rmsd_ = tag->getOption<Real>("max_overlap_rmsd", 4.0);
	attach_pose_ = core::pose::PoseOP( new pose::Pose );
	std::string fname( tag->getOption< std::string >( "attach_pdb" ) );


	core::import_pose::pose_from_file(*attach_pose_, fname , core::import_pose::PDB_file);
	minimize_after_overlap_ = tag->getOption<bool>("minimize_after_overlap",true);
	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}
	sfxn_->set_weight( scoring::linear_chainbreak, 2.0 );
}

std::string
MergePDBatOverlapMover::get_name() const {
	return MergePDBatOverlapMover::mover_name();
}

std::string
MergePDBatOverlapMover::mover_name()
{
	return "MergePDBatOverlapMover";
}

void MergePDBatOverlapMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction attachment_termini_type;
	attachment_termini_type.name( "attachment_termini_type" );
	attachment_termini_type.base_type( xs_string );
	attachment_termini_type.add_restriction( xsr_enumeration, "n_term" );
	attachment_termini_type.add_restriction( xsr_enumeration, "c_term" );
	xsd.add_top_level_element( attachment_termini_type );

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "attachment_termini", "attachment_termini_type", "XRW TO DO", "n_term" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_length", xsct_non_negative_integer, "XRW TO DO", "40" )
		+ XMLSchemaAttribute::attribute_w_default( "max_overlap_rmsd", xsct_real, "XRW TO DO", "2.0" )
		+ XMLSchemaAttribute::required_attribute( "attach_pdb", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "do_minimize", xsct_rosetta_bool, "Perform energy minimization", "true" )
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function used for packing and design." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Makes junctions from file", attlist );

}

} // simple_moves
} // protocols

