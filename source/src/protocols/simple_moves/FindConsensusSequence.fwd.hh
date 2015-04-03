// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/FindConsensusSequence.fwd.hh
/// @brief Finds the consensus sequence that optimizes energy of all input poses
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_FindConsensusSequence_fwd_hh
#define INCLUDED_protocols_simple_moves_FindConsensusSequence_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace simple_moves {

// Forward declaration of class
class FindConsensusSequence;

typedef utility::pointer::shared_ptr< FindConsensusSequence > FindConsensusSequenceOP;
typedef utility::pointer::shared_ptr< FindConsensusSequence const > FindConsensusSequenceCOP;

} // namespace simple_moves
} // namespace protocols

#endif //INCLUDED_protocols_simple_moves_FindConsensusSequence_fwd_hh
