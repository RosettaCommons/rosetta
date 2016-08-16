// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/medal/util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MEDAL_UTIL_HH
#define INCLUDED_PROTOCOLS_MEDAL_UTIL_HH

// External headers
#include <boost/unordered_set.hpp>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

namespace protocols {
namespace medal {

typedef utility::vector1<double> Probabilities;

/// @brief Computes per-residue alignment probabilities. To minimize assumptions
/// about the fold tree, the number of residues in the pose is an explicit input.
void alignment_probabilities(const unsigned num_residues,
	const core::sequence::SequenceAlignment& alignment,
	Probabilities* p);

/// @brief Computes per-residue chunk probabilities. A linear voting scheme is
/// used to bias sampling toward larger chunks.
void chunk_probabilities(const protocols::loops::Loops& chunks, Probabilities* p);

/// @brief Computes per-residue cutpoint probabilities. Sampling is biased
/// toward residues proximal to a cutpoint.
void cutpoint_probabilities(const unsigned num_residues,
	const core::kinematics::FoldTree& tree,
	Probabilities* p);

/// @brief Biases sampling away from termini
void end_bias_probabilities(const unsigned num_residues, Probabilities* p);

/// @brief Zeroes out the probability of residues for which fragment insertion
/// would cross cutpoint boundaries.
void invalidate_residues_spanning_cuts(const core::kinematics::FoldTree& tree,
	const core::Size fragment_len,
	Probabilities* probs);

/// @brief Populates a hash table with the residues contained in loops
void as_set(protocols::loops::LoopsCOP loops, boost::unordered_set<core::Size>* s);

void to_centroid(core::pose::Pose* pose);

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_UTIL_HH_
