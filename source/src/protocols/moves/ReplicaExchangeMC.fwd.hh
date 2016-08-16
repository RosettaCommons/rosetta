// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/ReplicaExchangeMC.fwd.hh
/// @brief  implementing a Recplica Exchange Monte Carlo Mover
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_protocols_moves_ReplicaExchangeMC_fwd_hh
#define INCLUDED_protocols_moves_ReplicaExchangeMC_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers
namespace protocols {
namespace moves {

class ReplicaExchangeMC;
typedef utility::pointer::shared_ptr< ReplicaExchangeMC > ReplicaExchangeMC_OP;
typedef utility::pointer::shared_ptr< ReplicaExchangeMC const > ReplicaExchangeMC_COP;

} // moves
} // prot

#endif

