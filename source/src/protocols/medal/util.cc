// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/medal/util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/medal/util.hh>

// C/C++ headers
#include <cmath>
#include <set>

// External headers
#include <boost/unordered/unordered_set.hpp>

// Utility headers
#include <numeric/prob_util.hh>
#include <utility/iter_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

namespace protocols {
namespace medal {

typedef utility::vector1<double> Probabilities;
using core::Size;

/// @brief Lower sampling probability near termini
void end_bias_probabilities(const Size num_residues, Probabilities* p) {
	debug_assert(p);

	// Penalize 10% lowest, highest residues
	const Size offset_left = static_cast<Size>(std::ceil(0.1 * num_residues));
	const Size offset_right = num_residues - offset_left + 1;

	p->clear();
	for ( Size i = 1; i <= num_residues; ++i ) {
		bool near_terminus = i <= offset_left || i >= offset_right;
		p->push_back(near_terminus ? 0.2 : 0.8);
	}
	numeric::normalize(p->begin(), p->end());
}

void alignment_probabilities(const Size num_residues,
	const core::sequence::SequenceAlignment& alignment,
	Probabilities* p) {
	debug_assert(p);

	p->clear();
	core::id::SequenceMapping mapping = alignment.sequence_mapping(1, 2);
	for ( Size i = 1; i <= num_residues; ++i ) {
		p->push_back(mapping[i] > 0 ? 0.2 : 0.8);
	}
	numeric::normalize(p->begin(), p->end());
}

/// @detail Linear voting based on chunk length
/// Assumes <chunks> are in ascending order
void chunk_probabilities(const protocols::loops::Loops& chunks, Probabilities* p) {
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	debug_assert(p);

	p->clear();
	for ( auto const & chunk : chunks ) {
		for ( Size j = chunk.start(); j <= chunk.stop(); ++j ) {
			p->push_back(chunk.length());
		}
	}
	numeric::normalize(p->begin(), p->end());
}

void cutpoint_probabilities(const Size num_residues, const core::kinematics::FoldTree& tree, Probabilities* p) {
	debug_assert(p);

	// List of cutpoints sorted by position
	std::set<Size> cutpoints;
	for ( Size i = 1; i <= tree.num_cutpoint(); ++i ) {
		Size cutpoint = tree.cutpoint(i);
		if ( cutpoint == num_residues ) {  // cutpoint at end of chain
			continue;
		}

		cutpoints.insert((Size)cutpoint);
	}

	p->clear();
	for ( Size residue = 1; residue <= num_residues; ++residue ) {
		const double nearest_cutpoint = *utility::find_closest(cutpoints.begin(), cutpoints.end(), residue);
		const double distance = std::abs(nearest_cutpoint - residue);
		p->push_back(std::pow(distance + 1, -3));
	}
	numeric::normalize(p->begin(), p->end());
}

/// @detail Given a fold tree and fragment library, zero out residues whose modification
/// would span a cutpoint boundary
void invalidate_residues_spanning_cuts(const core::kinematics::FoldTree& tree,
	const core::Size fragment_len,
	Probabilities* probs) {
	for ( Size i = 1; i <= tree.num_cutpoint(); ++i ) {
		const Size cutpoint = tree.cutpoint(i);
		for ( Size j = (cutpoint - fragment_len + 2); j <= cutpoint; ++j ) {
			(*probs)[j] = 0;
		}
	}
}

void as_set(protocols::loops::LoopsCOP loops, boost::unordered_set<core::Size>* s) {
	debug_assert(s);
	debug_assert(loops);

	s->clear();

	for ( auto const & i : *loops ) {
		for ( core::Size j = i.start(); j <= i.stop(); ++j ) {
			s->insert(j);
		}
	}
}

void to_centroid(core::pose::Pose* pose) {
	debug_assert(pose);
	if ( !pose->is_centroid() ) {
		core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID_t);
	}
}

}  // namespace medal
}  // namespace protocols
