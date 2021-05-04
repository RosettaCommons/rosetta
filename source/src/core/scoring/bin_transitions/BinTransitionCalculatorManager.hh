// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/bin_transitions/BinTransitionCalculatorManager
/// @brief A static singleton that ensures that BinTransitionCalculators load data from disk once,
/// lazily, and in a threadsafe manner.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_core_scoring_bin_transitions_BinTransitionCalculatorManager_hh
#define INCLUDED_core_scoring_bin_transitions_BinTransitionCalculatorManager_hh

// Unit headers
#include <core/scoring/bin_transitions/BinTransitionCalculatorManager.fwd.hh>
#include <core/scoring/bin_transitions/BinTransitionCalculator.fwd.hh>

// Utility header
#include <utility/SingletonBase.hh>
#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif //MULTI_THREADED

// C++ header
#include <map>
#include <string>

namespace core {
namespace scoring {
namespace bin_transitions {

/// @brief A static singleton that ensures that BinTransitionCalculators load data from disk once,
/// lazily, and in a threadsafe manner.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class BinTransitionCalculatorManager : public utility::SingletonBase< BinTransitionCalculatorManager > {
	friend class utility::SingletonBase< BinTransitionCalculatorManager >;

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	BinTransitionCalculatorManager();
	BinTransitionCalculatorManager(BinTransitionCalculatorManager const & ) = delete;
	BinTransitionCalculatorManager operator=(BinTransitionCalculatorManager const & ) = delete;

	/// @brief Create a BinTransitionCalculator from a file read from disk.  Used for onetime object creation; called
	/// ONLY by the BinTransationCalculatorManager and ONLY for one read of each unique file.
	static BinTransitionCalculatorOP create_bin_transition_calculator( std::string const & bin_transitions_file );

public: // Public methods //////////////////////////////////////////////////

	/// @brief Given a bin transition file, get a bin transition calculator.
	/// @details If we've already loaded the file and created a calculator from it, then we have cached the
	/// calculator and can return a clone to it.  Otherwise, we create the calculator, cache it, and return a
	/// clone to it.
	/// @note Fully threadsafe.
	BinTransitionCalculatorCOP get_bin_transition_calculator( std::string const & bin_transitions_file ) const;

private:  // Private data /////////////////////////////////////////////////////

#ifdef MULTI_THREADED
	mutable utility::thread::ReadWriteMutex mutex_;
#endif //MULTI_THREADED

	/// @brief All the BinTransitionCalculator objects that we've cached, with filenames as keys.
	mutable std::map< std::string, BinTransitionCalculatorOP > cached_bin_transition_calculators_;

};

} //bin_transitions
} //scoring
} //core

#endif //INCLUDED_core/scoring/bin_transitions_BinTransitionCalculatorManager_fwd_hh



