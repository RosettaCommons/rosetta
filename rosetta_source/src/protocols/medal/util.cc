// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/medal/util.hh>

// C/C++ headers
#include <cmath>
#include <set>

// Utility headers
#include <numeric/prob_util.hh>
#include <utility/iter_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

namespace protocols {
namespace medal {

typedef utility::vector1<double> Probabilities;

void alignment_probabilities(const unsigned num_residues,
														 const core::sequence::SequenceAlignment& alignment,
														 Probabilities* p) {
  assert(p);

  p->clear();
  core::id::SequenceMapping mapping = alignment.sequence_mapping(1, 2);
  for (unsigned i = 1; i <= num_residues; ++i) {
    p->push_back(mapping[i] > 0 ? 0.2 : 0.8);
  }
  numeric::normalize(p->begin(), p->end());
}

/// @detail Linear voting based on chunk length
/// Assumes <chunks> are in ascending order
void chunk_probabilities(const protocols::loops::Loops& chunks, Probabilities* p) {
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  assert(p);

  p->clear();
  for (Loops::const_iterator i = chunks.begin(); i != chunks.end(); ++i) {
    const Loop& chunk = *i;
    for (unsigned j = chunk.start(); j <= chunk.stop(); ++j) {
      p->push_back(chunk.length());
    }
  }
  numeric::normalize(p->begin(), p->end());
}

void cutpoint_probabilities(const unsigned num_residues, const core::kinematics::FoldTree& tree, Probabilities* p) {
  assert(p);

  // List of cutpoints sorted by position
	std::set<unsigned> cutpoints;
  for (unsigned i = 1; i <= tree.num_cutpoint(); ++i) {
    unsigned cutpoint = tree.cutpoint(i);
    if (cutpoint == num_residues)  // cutpoint at end of chain
      continue;

    cutpoints.insert(cutpoint);
  }

  p->clear();
  for (unsigned residue = 1; residue <= num_residues; ++residue) {
    const double nearest_cutpoint = *utility::find_closest(cutpoints.begin(), cutpoints.end(), residue);
    const double distance = std::abs(nearest_cutpoint - residue);
    p->push_back(std::pow(distance + 1, -1.5));
  }
  numeric::normalize(p->begin(), p->end());
}

/// @detail Given a fold tree and fragment library, zero out residues whose modification
/// would span a cutpoint boundary
void invalidate_residues_spanning_cuts(const core::kinematics::FoldTree& tree,
                                       const core::Size fragment_len,
                                       Probabilities* probs) {
  for (unsigned i = 1; i <= tree.num_cutpoint(); ++i) {
    const unsigned cutpoint = tree.cutpoint(i);
    for (unsigned j = (cutpoint - fragment_len + 2); j <= cutpoint; ++j) {
      (*probs)[j] = 0;
    }
  }
}

}  // namespace medal
}  // namespace protocols
