// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/HarmonicFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/func/HarmonicFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


#include <sstream>


// C++ Headers


namespace core {
namespace scoring {
namespace func {

FuncOP
HarmonicFunc::clone() const
{
	return FuncOP( new HarmonicFunc( x0_, sd_ ) );
}

Real
HarmonicFunc::func( Real const x ) const
{
	Real const z = ( x-x0_ )/sd_;
	return z * z;
}

Real
HarmonicFunc::dfunc( Real const x ) const
{
	return 2 * (x-x0_) / ( sd_ * sd_ );
}

void
HarmonicFunc::read_data( std::istream& in ) {
	in >> x0_ >> sd_;
}

void
HarmonicFunc::show_definition( std::ostream &out ) const {
	out << "HARMONIC " << x0_ << " " << sd_ << std::endl;
}

Size
HarmonicFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if (verbose_level > 100 ) {
		out << "HARM " <<  ( x - x0_ )/sd_ << std::endl;
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
