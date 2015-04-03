// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SigmoidFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Hetu Kamisetty


#include <core/scoring/func/SigmoidFunc.hh>
#include <math.h>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


#include <sstream>


// C++ Headers


namespace core {
namespace scoring {
namespace func {

Real
SigmoidFunc::func( Real const x ) const
{
	Real const z = slope_*( x-x0_ );
	return (1/(1+exp(-z)) - 0.5); // so at x0_, you get 0
}

Real
SigmoidFunc::dfunc( Real const x ) const
{
	Real fval = 1/(1+exp(-slope_*( x-x0_ )));
	return slope_*fval*fval*exp(-slope_*(x-x0_));
}

void
SigmoidFunc::read_data( std::istream& in ) {
	in >> x0_ >> slope_;
}

void
SigmoidFunc::show_definition( std::ostream &out ) const {
	out << "SIGMOID " << x0_ << " " << slope_ << std::endl;
}

Size
SigmoidFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	//I don't know what to put here.
	if (verbose_level > 100 ) {
		out << "HARM " <<  ( x - x0_ )/slope_ << std::endl;
	} else if (verbose_level > 70 ) {
		if ( x < x0_  && ( this->func(x) > threshold ) ) out << "-";
		else if ( x > x0_ && ( this->func(x) > threshold )) out << "+";
		else out << ".";
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core
