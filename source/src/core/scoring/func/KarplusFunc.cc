// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/KarplusFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author Nikolas Sgourakis

#include <core/scoring/func/KarplusFunc.hh>
#include <core/types.hh>

#include <iostream>

#include <cmath>


namespace core {
namespace scoring {
namespace func {

Real
KarplusFunc::func( Real const x ) const {
	Real const cosine = cos( x + Dphi_ );
	Real const j = A_ * cosine * cosine + B_ * cosine + C_;
	Real const z = (j - x0_)  / sd_;
	return z * z;

}

Real
KarplusFunc::dfunc( Real const x ) const {
	Real const sine = sin (x+Dphi_);
	Real const cosine = cos (x+Dphi_);
	Real const j = A_ * cosine * cosine + B_ * cosine + C_;
	return  -2 *( (j -x0_ ) / (sd_ * sd_) ) * sine *(2*A_*cosine + B_ ) ;
}

void
KarplusFunc::read_data( std::istream & in ) {
	in >> A_ >> B_ >>C_ >> Dphi_ >> x0_ >> sd_ >>offset_;
}

void KarplusFunc::show_definition( std::ostream & out ) const
{
	out << "KarplusFunc " << A_ <<' ' <<  B_ << ' ' << C_ << ' ' << Dphi_ << ' '<< x0_ << ' ' << sd_ << ' ' << offset_;
}

} // namespace constraints
} // namespace scoring
} // namespace core
