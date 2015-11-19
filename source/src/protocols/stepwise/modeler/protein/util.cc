// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinUtil
/// @brief a few functions used by several StepWiseProteinAnsatz classes
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/util.hh>

//////////////////////////////////
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/types.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>


#include <utility/numbers.hh>

using core::Real;
using core::Size;
using core::pose::Pose;
using utility::tools::make_vector1;
using ObjexxFCL::string_of;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

//////////////////////////////////////////////////////////////////////////////////////////////////
// void
// setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
// 	core::pose::Pose const & pose_2,
// 	Size const base_res,
// 	core::id::AtomID_Map< core::id::AtomID > & atom_ID_map){
// 	setup_protein_backbone_atom_id_map( pose_1, pose_2, base_res, atom_ID_map );
// }

//////////////////////////////////////////////////////////////////////////////////////////////////
void
setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
	core::pose::Pose const & pose_2,
	Size const base_res,
	Size const base_res2,
	core::id::AtomID_Map< core::id::AtomID > & atom_ID_map){

	using namespace core::id;

	if ( base_res == 0 ) return;
	if ( base_res2 == 0 ) return;

	if ( !pose_1.residue_type( base_res ).is_protein() ) return;
	if ( !pose_2.residue_type( base_res2 ).is_protein() ) return;

	if ( pose_1.residue_type( base_res ).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) return;
	if ( pose_2.residue_type( base_res2 ).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) return;

	{
		AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res ), pose_1 ));
		AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res2 ), pose_2 ));
		atom_ID_map.set( atom1, atom2 );
	}

	{
		AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", base_res ), pose_1 ));
		AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", base_res2 ), pose_2 ));
		atom_ID_map.set( atom1, atom2 );
	}

	{
		AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", base_res ), pose_1 ));
		AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", base_res2 ), pose_2 ));
		atom_ID_map.set( atom1, atom2 );
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_protein_modeling_info( pose::Pose const & pose,
	Size const moving_res,
	utility::vector1< Size > & moving_res_list ){

	if ( !pose.residue( moving_res ).is_protein() ) return;

	// go back another residue -- this was the default choice in protein SWA.
	utility::vector1< Size > const sample_res_for_pose = core::pose::full_model_info::get_sample_res_for_pose( pose );
	Size const upstream_res = pose.fold_tree().get_parent_residue( moving_res );
	if ( pose.residue_type( upstream_res ).is_protein() &&
			sample_res_for_pose.has_value( upstream_res ) ) moving_res_list.push_back( upstream_res );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_bridge_res( pose::Pose const & pose,
	utility::vector1< Size > const & moving_res_list /*working*/ ){

	// check if CCD-closure  is necessary.
	utility::vector1< Size > bridge_res;
	if ( moving_res_list.size() == 0 ) return bridge_res;

	Size const moving_res = moving_res_list[ 1 ];
	utility::vector1< Size > const cutpoints_closed = figure_out_moving_cutpoints_closed_from_moving_res( pose, moving_res );
	if ( cutpoints_closed.size() == 0 ) return bridge_res;

	for ( Size n = 1; n <= cutpoints_closed.size(); n++ ) {

		int cutpoint_closed = static_cast<int>( cutpoints_closed[ n ] );
		if ( !pose.residue_type( cutpoint_closed ).is_protein() ) continue;

		// as in protein SWA, choose two bridge residues for CCD closure on the 'other side' of sampled residue.
		utility::vector1< int > offsets;
		utility::vector1< Size > const sample_res_for_pose = core::pose::full_model_info::get_sample_res_for_pose( pose );
		if ( moving_res_list.has_value( cutpoint_closed ) ) {
			offsets = make_vector1( +1, +2 );
		} else if ( moving_res_list.has_value( cutpoint_closed+1 ) ) {
			offsets = make_vector1( -1, 0 );
		} else {
			offsets = make_vector1( 0, +1 ); // bracket CCD closure point, since moving_res does not.
		}
		utility::vector1< Size > working_bridge_res;
		for ( Size n = 1; n <= offsets.size(); n++ ) {
			int const bridge_res = cutpoint_closed + offsets[n];
			if ( bridge_res < 1 ) continue;
			runtime_assert( !moving_res_list.has_value( bridge_res ) );
			if ( sample_res_for_pose.has_value( bridge_res ) ) working_bridge_res.push_back( bridge_res );
		}
		bridge_res = merge_vectors( bridge_res, const_full_model_info( pose ).sub_to_full( working_bridge_res ) );
	}
	return bridge_res;

}


//////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
just_protein( utility::vector1< Size > const & res_list, pose::Pose const & pose ){
	utility::vector1< Size > protein_res_list;
	for ( Size n = 1; n <= res_list.size(); n++ ) { if ( pose.residue_type( res_list[n] ).is_protein() ) protein_res_list.push_back( res_list[n] ); }
	return protein_res_list;
}

//////////////////////////////////////////////////////////////////////////
bool
contains_protein( core::pose::Pose const & pose ){
	utility::vector1< Size > protein_res_list;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) { if ( pose.residue_type( n ).is_protein() ) return true; }
	return false;
}


} //protein
} //modeler
} //stepwise
} //protocols
