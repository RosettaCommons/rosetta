// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/assert.hh
/// @author Sergey Lyskov
///
/// @note Some assert's related macros


#ifndef INCLUDED_utility_assert_hh
#define INCLUDED_utility_assert_hh

#include <utility/backtrace.hh> // for debug_assert()

/// @brief Macro wrapper for paramters that used only indebug_assert(...) statments.
///        Intended to supress 'unused parameter' warning.
///        Ported from Rosetta++::Pack.cc
///
/// Example of usage: ResfileReader::read_aa_list(utility::vector1< std::string > const & ASSERT_ONLY(tokens) ) { ...}
///
#ifndef NDEBUG // Debug version
#define ASSERT_ONLY(x) x
#else
	#define ASSERT_ONLY(x)
#endif

#ifndef USEMPI
#define MPI_ONLY(x)
#else
  #define MPI_ONLY(x) x
#endif

namespace utility {
} // namespace utility

#endif // INCLUDED_utility_assert_HH
