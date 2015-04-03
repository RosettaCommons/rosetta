// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/FlatHarmonicFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Chris King modified from James Thompson


#include <core/scoring/func/FlatHarmonicFunc.hh>

#include <core/types.hh>
#include <math.h>

#include <utility/pointer/ReferenceCount.hh>


#include <sstream>


// C++ Headers


namespace core {
namespace scoring {
namespace func {

Real
FlatHarmonicFunc::func( Real const x ) const
{
		if( fabs( x - x0_ ) <= tol_ ) return 0.0;
		else if( x > x0_ ){
				Real const z = ( x - x0_ - tol_ ) / sd_;
				return z * z;
		}
		else{
				Real const z = ( x - x0_ + tol_ ) / sd_;
				return z * z;
		}
}

Real
FlatHarmonicFunc::dfunc( Real const x ) const
{
		if( fabs( x - x0_ ) <= tol_ ) return 0.0;
		else if( x > x0_ ) return 2 * ( x - x0_ - tol_ ) / ( sd_ * sd_ );
		else return 2 * ( x - x0_ + tol_ ) / ( sd_ * sd_ );
}

void
FlatHarmonicFunc::read_data( std::istream& in ) {
	in >> x0_ >> sd_ >> tol_;
}

void
FlatHarmonicFunc::show_definition( std::ostream &out ) const {
	out << "FLAT_HARMONIC " << x0_ << " " << sd_ << " " << tol_ << std::endl;
}

Size
FlatHarmonicFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if (verbose_level > 100 ) {
		if( fabs( x - x0_ ) <= tol_ ) out << "HARM " << 0.0 << std::endl;
		else if( x > x0_ ){
				Real const z = ( x - x0_ - tol_ ) / sd_;
				out << "FLHARM " << z << std::endl;
		}
		else{
				Real const z = ( x - x0_ + tol_ ) / sd_;
				out << "FLHARM " << z << std::endl;
		}
	} else if (verbose_level > 70 ) {
		if ( x < x0_ + tol_  && ( this->func(x) > threshold ) ) out << "-";
		else if ( x > x0_ + tol_ && ( this->func(x) > threshold )) out << "+";
		else out << ".";
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core
