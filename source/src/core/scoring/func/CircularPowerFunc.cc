// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularPowerFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <core/scoring/func/CircularPowerFunc.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

CircularPowerFunc::CircularPowerFunc(
	Real const x0_radians,
	Real const sd_radians,
	int const power,
	Real const weight
) :
	x0_( x0_radians ),
	sd_( sd_radians ),
	power_( power ),
	weight_( weight )
{
	assert( std::abs( std::pow(3.0,2) - 9.0 ) < 1e-3 );
	assert( power_ != 1 && power_ != 0 );
}

FuncOP
CircularPowerFunc::clone() const {
	return FuncOP( new CircularPowerFunc( *this ) );
}

Real
CircularPowerFunc::func( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return weight_ * std::pow( z, power_ );
}

Real
CircularPowerFunc::dfunc( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return weight_ * power_ * std::pow( z, power_ - 1 ) / sd_;
}

} // constraints
} // scoring
} // core
