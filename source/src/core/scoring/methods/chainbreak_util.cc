// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/chainbreak_util.cc
/// @brief  Utility functions for scoring chainbreaks.
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <boost/unordered/unordered_set.hpp>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

bool is_lower_cutpoint(
	core::Size residue,
	core::pose::Pose const & pose
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	const bool is_cutpoint_in_tree_lower  = pose.fold_tree().is_cutpoint(residue);
	const bool use_pose_cutpoint_variants = option[OptionKeys::score::score_pose_cutpoint_variants]();
	const bool has_lower_variant_type     = pose.residue(residue).has_variant_type(chemical::CUTPOINT_LOWER);
	return (has_lower_variant_type && (is_cutpoint_in_tree_lower || use_pose_cutpoint_variants));
}

bool is_upper_cutpoint(
	core::Size residue,
	core::pose::Pose const & pose
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	const bool is_cutpoint_in_tree_upper  = pose.fold_tree().is_cutpoint(residue - 1);
	const bool use_pose_cutpoint_variants = option[OptionKeys::score::score_pose_cutpoint_variants]();
	const bool has_upper_variant_type     = pose.residue(residue).has_variant_type(chemical::CUTPOINT_UPPER);
	if ( residue <= 1 ) return false;
	return (has_upper_variant_type && (is_cutpoint_in_tree_upper || use_pose_cutpoint_variants));
}

void find_cutpoint_variants(
	const core::pose::Pose& pose,
	const core::kinematics::FoldTree&,
	utility::vector1<int>* cutpoints
) {
	using core::Size;
	using boost::unordered_set;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	debug_assert(cutpoints);
	unordered_set<int> unique_cutpoints;

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( is_lower_cutpoint(ii,pose) ) unique_cutpoints.insert(ii);
	}

	// Update output parameter
	std::copy(unique_cutpoints.begin(),
		unique_cutpoints.end(),
		std::back_inserter(*cutpoints));

	std::sort(cutpoints->begin(), cutpoints->end());
}

} // namespace methods
} // namespace scoring
} // namespace core
