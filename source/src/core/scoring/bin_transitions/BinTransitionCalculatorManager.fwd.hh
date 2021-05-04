// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/bin_transitions/BinTransitionCalculatorManager.fwd.hh
/// @brief A static singleton that ensures that BinTransitionCalculators load data from disk once,
/// lazily, and in a threadsafe manner.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_core_scoring_bin_transitions_BinTransitionCalculatorManager_fwd_hh
#define INCLUDED_core_scoring_bin_transitions_BinTransitionCalculatorManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace scoring {
namespace bin_transitions {

class BinTransitionCalculatorManager;

using BinTransitionCalculatorManagerOP = utility::pointer::shared_ptr< BinTransitionCalculatorManager >;
using BinTransitionCalculatorManagerCOP = utility::pointer::shared_ptr< BinTransitionCalculatorManager const >;

} //bin_transitions
} //scoring
} //core

#endif //INCLUDED_core_scoring_bin_transitions_BinTransitionCalculatorManager_fwd_hh
