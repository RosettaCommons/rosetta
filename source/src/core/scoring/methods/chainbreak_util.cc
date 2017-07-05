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
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.chainbreak_util" );

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
	const bool is_cutpoint_in_tree_upper  = residue == 1 || pose.fold_tree().is_cutpoint(residue - 1);
	const bool use_pose_cutpoint_variants = option[OptionKeys::score::score_pose_cutpoint_variants]();
	const bool has_upper_variant_type     = pose.residue(residue).has_variant_type(chemical::CUTPOINT_UPPER);
	return (has_upper_variant_type && (is_cutpoint_in_tree_upper || use_pose_cutpoint_variants));
}

void find_cutpoint_variants(
	core::pose::Pose const & pose,
	core::kinematics::FoldTree const &,
	utility::vector1<int> & cutpoints
) {
	using core::Size;
	using boost::unordered_set;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	unordered_set<int> unique_cutpoints;

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( is_lower_cutpoint(ii,pose) ) unique_cutpoints.insert(ii);
	}

	// Update output parameter
	std::copy(unique_cutpoints.begin(),
		unique_cutpoints.end(),
		std::back_inserter(cutpoints));

	std::sort(cutpoints.begin(), cutpoints.end());
}

/// helper function for looking at residue connections to get lower/upper partners
bool
lower_upper_connected_across_cutpoint( core::conformation::Residue const & lower_rsd,
	core::conformation::Residue const & upper_rsd )
{
	for ( Size k = 1; k <= lower_rsd.connect_map_size(); k++ ) {
		if ( lower_rsd.residue_connect_atom_index( k ) != lower_rsd.upper_connect_atom() ) continue;
		Size upper( lower_rsd.connected_residue_at_resconn( k ) );
		if ( upper != upper_rsd.seqpos() ) return false;
		Size const m = lower_rsd.residue_connection_conn_id( k );
		if ( upper_rsd.residue_connect_atom_index( m ) != upper_rsd.lower_connect_atom() ) return false;
		runtime_assert( upper_rsd.connected_residue_at_resconn( m ) == lower_rsd.seqpos() );
		return true;
	}
	return false;
}

/// @details Instead of assuming cutpoint partner is simply cutpoint+1, find which residue connects via lower/upper.
///          Important in handling cyclized nucleotides.
Size
get_upper_cutpoint_partner_for_lower( pose::Pose const & pose, Size const lower_res )
{
	using namespace core::conformation;
	Residue const & lower_rsd( pose.residue( lower_res ) );
	// generally upper_cutpoint_partner is just lower_res + 1
	Size upper_cutpoint_partner( 0 );
	for ( Size k = 1; k <= lower_rsd.connect_map_size(); k++ ) {
		Size other( lower_rsd.connected_residue_at_resconn( k ) );
		if ( other == 0 ) continue;
		if ( lower_upper_connected_across_cutpoint( lower_rsd, pose.residue( other ) ) ) {
			upper_cutpoint_partner = other; break;
		}
	}

	return upper_cutpoint_partner;
}

Size
get_lower_cutpoint_partner_for_upper( pose::Pose const & pose, Size const upper_res )
{
	using namespace core::conformation;
	Residue const & upper_rsd( pose.residue( upper_res ) );
	// generally lower_cutpoint_partner is just upper_res - 1
	Size lower_cutpoint_partner( 0 );
	for ( Size k = 1; k <= upper_rsd.connect_map_size(); k++ ) {
		Size other( upper_rsd.connected_residue_at_resconn( k ) );
		if ( other == 0 ) continue;
		if ( lower_upper_connected_across_cutpoint( pose.residue( other ), upper_rsd ) ) {
			lower_cutpoint_partner = other; break;
		}
	}

	return lower_cutpoint_partner;
}

} // namespace methods
} // namespace scoring
} // namespace core
