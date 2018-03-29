// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/calc_taskop_movers/ConsensusDesignMover.fwd.hh
/// @brief fwd.hh file for ConsensusDesignMover
/// @author Florian Richter (floric@u.washington.edu), april 2011

#ifndef INCLUDED_protocols_simple_moves_ConsensusDesignMover_fwd_hh
#define INCLUDED_protocols_simple_moves_ConsensusDesignMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace calc_taskop_movers {

class ConsensusDesignMover;
typedef utility::pointer::shared_ptr< ConsensusDesignMover > ConsensusDesignMoverOP;

} // moves
} // protocols


#endif
