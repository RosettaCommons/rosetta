// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/util.cc
/// @brief  Miscellaneous functions that pop up everywhere in epr_deer calculation
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <core/scoring/epr_deer/util.hh>

#include <core/types.hh>

#include <numeric/constants.hh>

#include <boost/math/special_functions/erf.hpp>

#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.util" );

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
) {
	Real const top = std::exp( -0.5 * pow( ( x - avg ) / stdev, 2 ) );
	Real const bottom = stdev * sqrt( 2 * numeric::constants::d::pi );
	return top / bottom;
}

/// @brief Calculate probit function given a certain value
/// @param x: Value to calculate
/// @return Probit function of x
Real
probit(
	Real const & x
) {
	if ( x < 0.0 || x > 1.0 ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Probit function accessing number outside its range (0-1)!" );
	}
	return sqrt( 2.0 ) * boost::math::erf_inv( 2 * x - 1 );
}

/// @brief Calculate logit function given certain value
/// @param x: Value to calculate
/// @param lo: Lower bound
/// @param hi: Upper bound
/// @return: Logit value
Real
logit(
	Real const & x,
	Real const & lo, // = 0.0
	Real const & hi // = 1.0
) {
	Real const mod_x = ( x - lo ) / ( hi - lo );
	if ( mod_x < 0.0 || mod_x > 1.0 ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Logit function accessing number outside its range (0-1)!" );
	}
	return log( mod_x / ( 1.0 - mod_x ) );
}

/// @brief  Safe version of log() that doesn't crash when zero is passed
/// @param  q: Value to take logarithm
/// @return natural logarithm of q, or arbitrarily large number of q=0.0
Real
ln(
	Real const & q
) {
	if ( std::isinf( std::abs( log( q ) ) ) ) {
		if ( q == 0.0 ) {
			return log( std::numeric_limits< float >::min() );
		} else {
			return log( std::numeric_limits< float >::max() );
		}
	} else {
		return log( q );
	}
}

/// @brief  Return the minimum and maximum bins from two distributions
/// @param  histr1: First distribution
/// @param  histr2: Second distribution
/// @return Tuple of minimum and maximum, in that order
std::tuple< Size, Size >
min_max(
	std::map< Size, Real > const & histr1,
	std::map< Size, Real > const & histr2
) {
	return std::make_tuple(
		std::min( histr1.begin()->first,  histr2.begin()->first  ),
		std::max( histr1.rbegin()->first, histr2.rbegin()->first ) );
}

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
) {
	return min + ( max - min ) / ( 1 + std::exp( -1 * x ) );
}

}
}
}
