// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SquareWell2Func.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das
/// @author Jianqing Xu

#include <sstream>

#include <basic/Tracer.hh>
#include <core/scoring/func/SquareWell2Func.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

static thread_local basic::Tracer TR( "core.scoring.constraints.SquareWall2Func" );

// C++ Headers


namespace core {
namespace scoring {
namespace func {

Real
SquareWell2Func::func( Real const x ) const
{
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ );
	// TR<<"x="<<x<<"   x0_="<<x0_<<"   |x-x0|="<<std::fabs(z) <<"  x_range_="<<x_range_<<std::endl;
	if ( std::fabs(z) > x_range_ ) {
		return well_depth_;
	}
	return 0.0;
}

Real
SquareWell2Func::dfunc( Real const /*x*/ ) const
{
	return 0.0; //This is bad news for the minimizer...
}

void
SquareWell2Func::read_data( std::istream& in ) {
	in >> x0_ >> x_range_ >>well_depth_;
	// check for optional DEGREES tag
	while ( in.peek() == ' ' || in.peek() == '\t' ) { in.get(); }
	if ( in.peek() == 'D' ) {
		std::string degrees;
		in >> degrees;
		if ( degrees == "DEGREES" ) {
			x0_      *= numeric::constants::r::deg2rad;
			x_range_ *= numeric::constants::r::deg2rad;
		}
	}
}

void
SquareWell2Func::show_definition( std::ostream &out ) const {
	out << "SQUARE_WELL_2 " << x0_ << " " << x_range_ << " " << well_depth_ << std::endl;
}

Size
SquareWell2Func::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "SQUARE_WELL_2 " <<  ( x < x0_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core
