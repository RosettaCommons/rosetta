// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/random.functions.hh
/// @brief  Random number generation functions
/// @author Ion Yannopoulos (ion@rosettacommons.org)
///
/// @note   When designing RNGs consider:
///         http://www.sgi.com/tech/stl/RandomNumberGenerator.html


#ifndef INCLUDED_utility_random_functions_hh
#define INCLUDED_utility_random_functions_hh


// C++ headers
#include <cstdlib>


namespace utility {


/// @brief Generic random number generator
/// @note  Default to std::rand in the absence of a better alternative
/// @note  On POSIX systems (e.g. Linux) use drand48 for doubles
template< typename T >
inline
T
random()
{
	using std::rand;

	return static_cast< T >( rand() ) / static_cast< T >( RAND_MAX );
}


/// @brief Generate random doubles.
/// @note  This should probably use C stdlib.h's drand48 if available.
/// @deprecated  Use random<double> instead.
inline
double
drand()
{
	using std::rand;

	return (double)rand() / (double)RAND_MAX;
}


} // namespace utility


#endif // INCLUDED_utility_random_functions_HH
