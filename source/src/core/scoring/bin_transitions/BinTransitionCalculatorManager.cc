// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/bin_transitions/BinTransitionCalculatorManager.cc
/// @brief A static singleton that ensures that BinTransitionCalculators load data from disk once, lazily, and in a threadsafe manner.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <core/scoring/bin_transitions/BinTransitionCalculatorManager.hh>

// Unit headers
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>

// Project header
#include <core/types.hh>

// Utility headers
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <functional>

// Construct tracer.
static basic::Tracer TR( "core.scoring.bin_transitions.BinTransitionCalculatorManager" );

namespace core {
namespace scoring {
namespace bin_transitions {

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
BinTransitionCalculatorManager::BinTransitionCalculatorManager()
{}

/// @brief Create a BinTransitionCalculator from a file read from disk.  Used for onetime object creation; called
/// ONLY by the BinTransationCalculatorManager and ONLY for one read of each unique file.
/*static*/
BinTransitionCalculatorOP
BinTransitionCalculatorManager::create_bin_transition_calculator(
	std::string const & bin_transitions_file
) {
	return BinTransitionCalculatorOP( new BinTransitionCalculator(bin_transitions_file) ); //Huh.  I can't seem to use make_shared here.
}

// Public methods /////////////////////////////////////////////////////////////

/// @brief Given a bin transition file, get a bin transition calculator.
/// @details If we've already loaded the file and created a calculator from it, then we have cached the
/// calculator and can return a clone to it.  Otherwise, we create the calculator, cache it, and return a
/// clone to it.
/// @note Fully threadsafe.
BinTransitionCalculatorCOP
BinTransitionCalculatorManager::get_bin_transition_calculator(
	std::string const & bin_transitions_file
) const {
	std::function< BinTransitionCalculatorOP () > creator( std::bind( &BinTransitionCalculatorManager::create_bin_transition_calculator, std::cref(bin_transitions_file) ) );
	return (
		utility::thread::safely_check_map_for_key_and_insert_if_absent(
		creator, SAFELY_PASS_MUTEX( mutex_ ), bin_transitions_file, cached_bin_transition_calculators_
		)
		)->clone();
}

} //bin_transitions
} //scoring
} //core
