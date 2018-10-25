// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/indexed_structure_store/utility.cc

#include <core/types.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/VariantType.hh>

#include <protocols/indexed_structure_store/utility.hh>

namespace protocols { namespace indexed_structure_store {

core::pose::PoseOP
append_pose_by_bond (
	core::pose::Pose & upstream_pose,
	core::pose::Pose & downstream_pose
) {
	return append_pose_with_overlap(
		upstream_pose,
		downstream_pose,
		0,
		DELETE_UPSTREAM
	);
}

core::pose::PoseOP
append_pose_with_overlap(
	core::pose::Pose & upstream_pose,
	core::pose::Pose & downstream_pose,
	core::Size overlap, OverlapDeleteMode delete_overlap
) {
	runtime_assert( upstream_pose.size() > overlap );
	runtime_assert( downstream_pose.size() > overlap );


	using namespace core::pose;
	PoseOP result(new Pose(upstream_pose));
	core::Size join_point = upstream_pose.size();
	core::Size org_jumps = upstream_pose.num_jump();

	result->append_pose_by_jump(downstream_pose, result->size());

	if ( overlap > 0 ) {
		if ( delete_overlap == DELETE_UPSTREAM ) {
			result->delete_residue_range_slow(join_point - overlap + 1, join_point);
			join_point -= overlap;
		} else {
			result->delete_residue_range_slow(join_point + 1, join_point + overlap);
		}
	}

	remove_variant_type_from_pose_residue( *result, core::chemical::UPPER_TERMINUS_VARIANT, join_point );
	remove_variant_type_from_pose_residue( *result, core::chemical::CUTPOINT_LOWER, join_point );
	remove_variant_type_from_pose_residue( *result, core::chemical::CUTPOINT_UPPER, join_point );
	remove_variant_type_from_pose_residue( *result, core::chemical::LOWER_TERMINUS_VARIANT, join_point + 1 );
	remove_variant_type_from_pose_residue( *result, core::chemical::CUTPOINT_LOWER, join_point + 1 );
	remove_variant_type_from_pose_residue( *result, core::chemical::CUTPOINT_UPPER, join_point + 1 );

	result->conformation().delete_chain_ending(join_point);

	result->conformation().update_polymeric_connection(join_point, true);
	core::kinematics::FoldTree ft( result->fold_tree() );
	ft.delete_jump_and_intervening_cutpoint( org_jumps + 1 );
	result->fold_tree( ft );

	return result;
}

} }
