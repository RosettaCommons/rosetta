// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/PeriodicFunc.cc
/// @brief Definition for periodic functions
/// @author Florian Richter, floric@u.washington.edu

#include <core/scoring/func/PeriodicFunc.hh>
#include <math.h>
#include <iostream>

namespace core {
namespace scoring {
namespace func {


Real
PeriodicFunc::func( Real const x ) const
{
  return ( k_ * (cos( n_periodic_ * ( x - x0_ ) ) ) ) + C_;
}

Real
PeriodicFunc::dfunc( Real const x ) const
{
  return -1. * k_ * n_periodic_ * sin( n_periodic_ * (x - x0_ ) );
}

void
PeriodicFunc::read_data( std::istream& in )
{
  in >> x0_ >> n_periodic_ >> k_ >> C_;
}

void
PeriodicFunc::show_definition(std::ostream &out ) const {
  out << "PERIODIC " << x0_ << " " << n_periodic_ << " " << k_ << " " << C_
		<< std::endl;
}

//copied from HarmonicFunc.cc
Size
PeriodicFunc::show_violations(
	std::ostream& out, Real x, Size verbose_level, Real threshold
) const {
  if (verbose_level > 100 ) {
    out << "PERIODIC " <<  func(x) << std::endl;
  }
  return Func::show_violations( out, x, verbose_level, threshold);

}


} //constraints
} //scoring
} //core
