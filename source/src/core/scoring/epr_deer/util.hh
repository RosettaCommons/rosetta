// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/util.hh
/// @brief  Miscellaneous functions that pop up everywhere in epr_deer calculation
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <map>
#include <tuple>

#include <core/types.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Compute amplitude at a given point of a normal distribution
/// @param  x: The point of interest
/// @param avg: Average value of distribution
/// @param stdev: Width of distribution
/// @return Amplitude
Real
gauss(
	Real const & x,
	Real const & avg,
	Real const & stdev
);

/// @brief Calculate probit function given a certain value
/// @param x: Value to calculate
/// @return Probit function of x
Real
probit(
	Real const & x
);

/// @brief  Safe version of log() that doesn't crash when zero is passed
/// @param  q: Value to take logarithm
/// @return natural logarithm of q, or arbitrarily large number of q=0.0
Real
ln(
	Real const & q
);

/// @brief  Return the minimum and maximum bins from two distributions
/// @param  histr1: First distribution
/// @param  histr2: Second distribution
/// @return Tuple of minimum and maximum, in that order
std::tuple< Size, Size >
min_max(
	std::map< Size, Real > const & histr1,
	std::map< Size, Real > const & histr2
);

/// @brief Sigmoid function implementation
/// @param  x: Number to convert
/// @param  min: Minimum value of output
/// @param  max: Maximum value of output
/// @return Y output on sigmoid function
/// @detail Note that neither the slope nor the y-intercept can be modified
Real
sigmoid(
	Real const & x,
	Real const & min,
	Real const & max
);

/// @brief Simple fxn to avoid repeating this code everywhere
/// @param map: Simulated DEER histogram to modify
/// @param key: Bin/X-value
/// @param val: Y-value to either add or append
template <typename T1, typename T2>
void
add_to_map(
	std::map< T1, T2 > & map,
	T1 const & key,
	T2 const & val
) {
	if ( map.find( key ) == map.end() ) {
		map[ key ] = val;
	} else {
		map[ key ] += val;
	}
}

/// @brief Simple fxn to avoid repeating this code everywhere
/// @param map: Simulated DEER histogram to modify
/// @param key: Bin/X-value
/// @param val: Y-value to either add or append
template <typename T1, typename T2>
void
append_to_map(
	std::map< T1, utility::vector1< T2 > > & map,
	T1 const & key,
	T2 const & val
) {
	if ( map.find( key ) == map.end() ) {
		map[ key ] = utility::vector1< T2 >{ val };
	} else {
		map[ key ].push_back( val );
	}
}

}
}
}
