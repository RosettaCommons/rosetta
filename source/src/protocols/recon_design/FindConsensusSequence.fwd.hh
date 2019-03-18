// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recon_design/FindConsensusSequence.fwd.hh
/// @brief Finds the consensus sequence that optimizes energy of all input poses.
/// Used in RECON msd application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_recon_design_FindConsensusSequence_fwd_hh
#define INCLUDED_protocols_recon_design_FindConsensusSequence_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace recon_design {

// Forward declaration of class
class FindConsensusSequence;

typedef utility::pointer::shared_ptr< FindConsensusSequence > FindConsensusSequenceOP;
typedef utility::pointer::shared_ptr< FindConsensusSequence const > FindConsensusSequenceCOP;

} // namespace recon_design
} // namespace protocols

#endif //INCLUDED_protocols_recon_design_FindConsensusSequence_fwd_hh
