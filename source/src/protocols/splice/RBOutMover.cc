// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RBOutMover.cc
/// @brief

// Unit headers
#include <protocols/splice/RBOutMover.hh>
#include <protocols/splice/RBOutMoverCreator.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.simple_moves.RBOutMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/conformation/Residue.hh>
#include <protocols/simple_moves/CutChainMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <fstream>
#include <algorithm>
#include <vector>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/splice/util.hh>
#include <algorithm>

namespace protocols {
namespace splice {


using namespace::protocols;

// XRW TEMP std::string
// XRW TEMP RBOutMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RBOutMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RBOutMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RBOutMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RBOutMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "RBOut";
// XRW TEMP }

RBOutMover::RBOutMover(): moves::Mover("RBOut"),
	template_pdb_fname_(""),
	jump_dbase_fname_(""),
	jump_from_foldtree_( false )
{
}

RBOutMover::~RBOutMover() = default;

bool compare(std::pair<core::Size,core::Size> p1, std::pair<core::Size, core::Size> p2) {return p1.first == p2.first;}

/// find disulfide pairs. Simply look for pairs of SG atoms on residues that are disulfide variants and within 3.0A
utility::vector1< std::pair< core::Size, core::Size > >
find_disulfs_in_range( core::pose::Pose const & pose, core::Size const start, core::Size const end ){
	utility::vector1< std::pair< core::Size, core::Size > > disulfs;
	disulfs.clear();

	runtime_assert( end <= pose.size() && end >= 1 && start <= end && start >= 1 );
	for ( core::Size resi = start; resi < end; ++resi ) {
		for ( core::Size resj = resi + 1 ; resj <= end; ++resj ) {
			if ( pose.residue( resi ).has_variant_type( core::chemical::DISULFIDE ) && pose.residue( resj ).has_variant_type( core::chemical::DISULFIDE ) && pose.residue( resi ).xyz( "SG" ).distance( pose.residue( resj ).xyz("SG")) <= 3.0 ) {
				disulfs.push_back( std::pair< core::Size, core::Size >( resi, resj )) ;
			}
		}
	}

	utility::vector1< std::pair<core::Size,core::Size> >::iterator it;
	it = std::unique(disulfs.begin(),disulfs.end(), compare);
	disulfs.resize(std::distance(disulfs.begin(),it) );

	return(disulfs);
}


core::kinematics::Jump
RBOutMover::get_disulf_jump( Pose & pose, core::pose::Pose & template_pose )
{
	core::Size const pose_vl_vh = find_vl_vh_cut(pose);
	core::Size const template_vl_vh = find_vl_vh_cut(template_pose);
	utility::vector1<core::Size> template_cys_pos,pose_cys_pos;
	find_disulfide_postions(template_pose,template_cys_pos);
	find_disulfide_postions(pose,pose_cys_pos);
	TR<<"Template cys positions: "<<template_cys_pos<<std::endl;
	TR<<"pose cys positions: "<<pose_cys_pos<<std::endl;
	TR<<"Template vl_Vh cut: "<<template_vl_vh<<std::endl;
	TR<<"pose vl_Vh cut: "<<pose_vl_vh<<std::endl;

	// I assume the order of the template is VL/VH. The jump will be from the heavy chain second disulfide to the light chain second disulfide

	core::Real const dist1( template_pose.residue( template_cys_pos[4] ).xyz( "CA" ).distance( pose.residue( pose_cys_pos[4] ).xyz( "CA" ) ) );
	//core::Real const dist2( template_pose.residue( template_cys_pos[4] ).xyz( "CA" ).distance( pose.residue( pose_cys_pos[2] ).xyz( "CA" ) ) );

	core::Size pose_vh_second_cys, pose_vl_second_cys, pose_vh_first_cys, pose_vl_first_cys= 0;
	if ( dist1<3 ) {
		pose_vh_second_cys = pose_cys_pos[4];
		pose_vh_first_cys = pose_cys_pos[3];
		pose_vl_second_cys = pose_cys_pos[2];
		pose_vl_first_cys = pose_cys_pos[1];
	} else {
		pose_vh_second_cys = pose_cys_pos[2];
		pose_vh_first_cys = pose_cys_pos[1];
		pose_vl_second_cys = pose_cys_pos[4];
		pose_vl_first_cys = pose_cys_pos[3];
	}
	if ( (pose_vh_second_cys==0)||(pose_vl_second_cys==0) ) {
		utility_exit_with_message("Did not find pose aligned cys. Check alignment between pose and template\n");
	}

	utility::vector1<std::array<int, 3>> template_ft_nodes(set_fold_tree_nodes(template_pose,template_cys_pos, template_vl_vh));

	core::kinematics::FoldTree template_ft;
	template_ft.clear();

	for ( auto i: template_ft_nodes ) {
		TR<<i[0]<<","<<i[1]<<","<<i[2]<<std::endl;
		template_ft.add_edge(i[0],i[1],i[2]);
	}
	template_ft.add_edge((int)template_cys_pos[4],(int)template_cys_pos[2],1);
	template_ft.reorder((int)template_cys_pos[3]);
	template_ft.delete_self_edges();
	template_ft.check_fold_tree();
	template_pose.fold_tree( template_ft );

	TR<<"pose fold tree:"<<pose.fold_tree()<<std::endl;
	TR<<"template_pose fold tree:"<<template_pose.fold_tree()<<std::endl;


	utility::vector1<core::Size> template_VH_cys_pos;
	utility::vector1<core::Size> template_VL_cys_pos;
	template_pose.conformation().insert_chain_ending(template_vl_vh);
	core::pose::PoseOP template_VH = template_pose.split_by_chain(2);
	core::pose::PoseOP template_VL = template_pose.split_by_chain(1);
	find_disulfide_postions(*template_VH,template_VH_cys_pos);
	find_disulfide_postions(*template_VL,template_VL_cys_pos);

	superimpose_source_on_pose( pose,pose_vh_first_cys,pose_vh_second_cys, *template_VH,template_VH_cys_pos[1],template_VH_cys_pos[2]);
	superimpose_source_on_pose( pose,pose_vl_first_cys,pose_vl_second_cys, *template_VL,template_VL_cys_pos[1],template_VL_cys_pos[2]);



	template_VL->append_pose_by_jump(*template_VH,template_VL_cys_pos[2]);
	template_VL->conformation().delete_chain_ending(template_VL->conformation().chain_end(1));

	template_VL->fold_tree(template_ft);
	TR<<"Fold Tree of template after alignmnet:"<<template_VL->fold_tree()<<std::endl;
	TR<<"Template jump before change:"<<template_pose.jump( 1 )<<std::endl;
	TR<<"Template jump after change:"<<template_VL->jump( 1 )<<std::endl;

	//template_VL->dump_pdb("template_after_combine_VL_VH.pdb");

	template_pose.set_jump( 1,template_VL->jump( 1 ) );
	//template_pose.dump_pdb("after_ft_change.pdb");
	return template_VL->jump( 1 );
}

void
RBOutMover::apply( Pose & pose )
{
	core::pose::Pose template_pose;
	core::import_pose::pose_from_file( template_pose, template_pdb_fname() , core::import_pose::PDB_file);

	core::kinematics::Jump const pose_disulf_jump = ( jump_from_foldtree() ? pose.jump( 1 ) : get_disulf_jump( pose, template_pose ) );
	//write to file
	std::string const source_pdb_name(parse_pdb_code(pose.pdb_info()->name()));
	std::ofstream data;
	data.open( jump_dbase_fname().c_str(), std::ios::out );
	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open jump dbase file for writing: " + jump_dbase_fname() + "\n" );
	}
	data<<pose_disulf_jump<<"\t"<<source_pdb_name<<std::endl;
	TR<<"Wrote jump info: "<<pose_disulf_jump<<source_pdb_name<<std::endl;
}

// XRW TEMP std::string
// XRW TEMP RBOutMover::get_name() const {
// XRW TEMP  return RBOutMover::mover_name();
// XRW TEMP }

moves::MoverOP
RBOutMover::clone() const
{
	return moves::MoverOP( new RBOutMover( *this ) );
}

moves::MoverOP
RBOutMover::fresh_instance() const
{
	return moves::MoverOP( new RBOutMover );
}

void
RBOutMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	template_pdb_fname( tag->getOption< std::string >( "template_fname" ) );
	jump_dbase_fname( tag->getOption< std::string >( "jump_dbase_fname" ));
	jump_from_foldtree( tag->getOption< bool >( "jump_from_foldtree", false ) );
	TR<<"Template pdb fname: "<<template_pdb_fname()<<" jump_dbase_fname: "<<jump_dbase_fname()<<" jump_from_foldtree: "<<jump_from_foldtree()<<std::endl;
}

std::string RBOutMover::get_name() const {
	return mover_name();
}

std::string RBOutMover::mover_name() {
	return "RBOut";
}

void RBOutMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "template_fname", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "jump_dbase_fname", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_from_foldtree", xsct_rosetta_bool, "XRW TO DO", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}



std::string RBOutMoverCreator::keyname() const {
	return RBOutMover::mover_name();
}

protocols::moves::MoverOP
RBOutMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RBOutMover );
}

void RBOutMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RBOutMover::provide_xml_schema( xsd );
}

void RBOutMover::find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos) {
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(core::chemical::DISULFIDE) ) {
			cys_pos.push_back(i);
		}
	}
}

utility::vector1<std::array<int, 3>> RBOutMover::set_fold_tree_nodes(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos, core::Size vl_vh_cut){
	utility::vector1<std::array<int, 3>> fold_tree_nodes;
	fold_tree_nodes.clear();

	/*Build the following Fold tree for the Vl/Vh fragment
	___________________________
	|                         |
	|                         |
	<------*<--------*------>//<-----*<-------*------->

	*/
	fold_tree_nodes.push_back({{(int) cys_pos[1],1,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[2],(int) cys_pos[1],-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[2],(int) vl_vh_cut,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[3],(int) vl_vh_cut+1,-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[4],(int) cys_pos[3],-1}});
	fold_tree_nodes.push_back({{(int) cys_pos[4],(int) pose.conformation().chain_end(1),-1}});
	// fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[2],1}});
	return fold_tree_nodes;
}

core::Size RBOutMover::find_vl_vh_cut(core::pose::Pose pose)  {
	protocols::simple_moves::CutChainMover ccm;
	return ccm.chain_cut(pose); //
}

void
RBOutMover::superimpose_source_on_pose( core::pose::Pose const & target_pose,core::Size target_from_res,core::Size target_to_res, core::pose::Pose & template_pose ,
	core::Size template_from_res,core::Size template_to_res){
	using namespace protocols::toolbox;
	utility::vector1< core::Size > target_positions, template_positions;
	target_positions.clear(); template_positions.clear();
	core::Size max_target_to_res  = std::min(target_pose.total_residue(),target_to_res+10);
	core::Size max_template_to_res  = std::min(template_pose.total_residue(),template_to_res+10);
	core::Size min_target_to_res  = std::max((int)target_from_res,(int)target_to_res-10);
	core::Size min_template_to_res  = std::max((int)template_from_res,(int)template_to_res-10);

	core::Size max_target_from_res  = std::min(target_to_res,target_from_res+10);
	core::Size max_template_from_res  = std::min(template_to_res,template_from_res+10);
	core::Size min_target_from_res  = std::max(1,(int)target_from_res-10);
	core::Size min_template_from_res  = std::max(1,(int)template_from_res-10);


	for ( core::Size i=min_template_from_res; i<=max_template_from_res; ++i ) {
		template_positions.push_back( i );
	}
	for ( core::Size i=min_template_to_res; i<=max_template_to_res; ++i ) {
		template_positions.push_back( i );
	}
	for ( core::Size i=min_target_from_res; i<=max_target_from_res; ++i ) {
		target_positions.push_back( i );
	}
	for ( core::Size i=min_target_to_res; i<=max_target_to_res; ++i ) {
		target_positions.push_back( i );
	}
	TR<<"Aligning template positions:"<<template_positions<<std::endl;
	TR<<"To target positions:"<<target_positions<<std::endl;

	utility::vector1< numeric::xyzVector< core::Real > > target_coor( coords( target_pose, target_positions ) ),
		template_coor( coords( template_pose, template_positions ));
	TR<<"template ref coords: "<<template_coor[1][1]<<std::endl;
	TR<<"source ref coords: "<<target_coor[1][1]<<std::endl;

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	superposition_transform( template_coor,target_coor ,rotation, to_init_center, to_fit_center );

	apply_superposition_transform( template_pose, rotation, to_init_center, to_fit_center );
	/// DEBUGGING
	//template_pose.dump_pdb(template_pose.sequence(1,1)+"superimposed_source_pose.pdb" );
}


} // splice
} // protocols

