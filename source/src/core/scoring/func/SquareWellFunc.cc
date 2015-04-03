// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SquareWellFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das


#include <core/scoring/func/SquareWellFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


#include <utility/vector1.hh>
#include <sstream>


// C++ Headers


namespace core {
namespace scoring {
namespace func {

Real
SquareWellFunc::func( Real const x ) const
{
	if ( x < x0_ ) {
		return well_depth_;
	}
	return 0.0;
}

Real
SquareWellFunc::dfunc( Real const /*x*/ ) const
{
	return 0.0; //This is bad news for the minimizer...
}

void
SquareWellFunc::read_data( std::istream& in ) {
	in >> x0_ >> well_depth_;
}

void
SquareWellFunc::show_definition( std::ostream &out ) const {
	out << "SQUARE_WELL " << x0_ << " " << well_depth_ << std::endl;
}

Size
SquareWellFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if (verbose_level > 100 ) {
		out << "SQUARE_WELL " <<  ( x < x0_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core
