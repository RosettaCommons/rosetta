// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularHarmonicFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <numeric/angle.functions.hh>
#include <core/types.hh>

#include <iostream>

namespace core {
namespace scoring {
namespace func {

Real
CircularHarmonicFunc::func( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return z * z + offset_;
}

Real
CircularHarmonicFunc::dfunc( Real const x ) const {
	return 2 * ( numeric::nearest_angle_radians(x,x0_)-x0_ ) / ( sd_ * sd_ );
}

void
CircularHarmonicFunc::read_data( std::istream & in ) {
	in >> x0_ >> sd_;
}

void CircularHarmonicFunc::show_definition( std::ostream & out ) const
{
	out << "CircularHarmonicFunc " << x0_ << ' ' << sd_;
}

} // namespace constraints
} // namespace scoring
} // namespace core
