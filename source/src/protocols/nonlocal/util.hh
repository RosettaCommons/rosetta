// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_UTIL_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_UTIL_HH

// C/C++ headers

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

// Package headers
#include <core/scoring/ScoreType.hh>

#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>


namespace protocols {
namespace nonlocal {

/// @brief Combine aligned and unaligned regions, limit size of final loop
protocols::loops::Loops combine_and_trim(core::Size min_chunk_sz,
	core::Size num_residues,
	const protocols::loops::LoopsOP aligned_regions,
	const protocols::loops::LoopsOP unaligned_regions);

// -- Utility methods -- not to be called directly
void find_regions_with_minimum_size(const core::sequence::SequenceAlignment& alignment,
	const core::Size unaligned_region_min_sz,
	protocols::loops::LoopsOP & aligned_regions,
	protocols::loops::LoopsOP & unaligned_regions);

/// @brief Best-effort attempt to limit the length of a chunk by recursively
/// decomposing <regions> such that min_chunk_sz <= |chunk| <= max_chunk_sz.
void limit_chunk_size(core::Size min_chunk_sz,
	core::Size max_chunk_sz,
	protocols::loops::LoopsOP & regions);

/// @brief Recursively decomposes <loop> into a series of <pieces>, each having
/// length less than or equal to <max_length>.
void decompose(core::Size min_chunk_sz,
	core::Size max_chunk_sz,
	const protocols::loops::Loop& loop,
	utility::vector1<protocols::loops::Loop>* pieces);

/// @brief If -abinitio:debug is enabled, writes <pose> to <file>.
void emit_intermediate(const core::pose::Pose& pose, const std::string& file);

/// @brief Computes the distance between consecutive CA atoms. If the distance exceeds
/// a user-specified threshold, creates a new chunk and adds it to <chunks>. CA-CA
/// distance threshold is retrieved from the option system (rigid::max_ca_ca_dist).
void chunks_by_CA_CA_distance(const core::pose::Pose& pose, protocols::loops::LoopsOP chunks);

/// @brief Computes the distance between consecutive CA atoms. If the distance exceeds
/// <threshold>, creates a new chunk and adds it to <chunks>.
void chunks_by_CA_CA_distance(const core::pose::Pose& pose, protocols::loops::LoopsOP chunks, double threshold);

/// @brief Returns the unweighted score of the ScoreType for the given residue. Assumes that the Pose has recently been scored by ScoreFunction with non-zero weight for the ScoreType.
core::Real get_per_residue_score(
	core::Size rsd_idx,
	core::scoring::ScoreType scoretype,
	core::pose::Pose const & pose
);

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_UTIL_HH_
