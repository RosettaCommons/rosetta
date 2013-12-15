// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/LinearPenaltyFunction.hh
/// @author Dominik Gront


#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/types.hh>

// C++ Headers
#include <cmath>
#include <iostream>


namespace core {
namespace scoring {
namespace func {

Real
LinearPenaltyFunction::func( Real const x ) const {

	Real dev = fabs(x-x_middle_);
	if ( dev <= half_width_ ) {
		return well_depth_;
	}
	return well_depth_ + slope_ * (dev-half_width_);
}

Real
LinearPenaltyFunction::dfunc( Real const /*x*/ ) const
{
	return 0.0;
}

void
LinearPenaltyFunction::read_data( std::istream& in ) {
	in >> x_middle_ >> well_depth_ >> half_width_ >> slope_;
}

void
LinearPenaltyFunction::show_definition( std::ostream &out ) const {
	out << "LINEAR_PENALTY " << x_middle_ << " " << well_depth_ <<  " " << half_width_ <<
	    " " << slope_ << std::endl;
}

Size
LinearPenaltyFunction::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {

	if (verbose_level > 100 ) {
		out << "LINEAR_PENALTY " <<  ( x < x_middle_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core
