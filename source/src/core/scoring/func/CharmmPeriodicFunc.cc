// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CharmmPeriodicFunc.cc
/// @brief Definition for periodic functions
/// @author Florian Richter, floric@u.washington.edu


#include <core/scoring/func/CharmmPeriodicFunc.hh>

#include <math.h>

#include <iostream>


namespace core {
namespace scoring {
namespace func {


Real
CharmmPeriodicFunc::func( Real const x ) const
{
  return 0.5 * k_ * (1 - cos( n_periodic_ * ( x - x0_ ) ) );
}

Real
CharmmPeriodicFunc::dfunc( Real const x ) const
{
  return 0.5 * k_ * n_periodic_ * sin( n_periodic_ * (x - x0_ ) );
}

void
CharmmPeriodicFunc::read_data( std::istream& in )
{
  in >> x0_ >> n_periodic_ >> k_;
}

void
CharmmPeriodicFunc::show_definition(std::ostream &out ) const
{
  out << "CHARMM_PERIODIC " << x0_ << " " << n_periodic_ << " " << k_ << std::endl;
}

//copied from HarmonicFunc.cc
Size
CharmmPeriodicFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const
{
  if (verbose_level > 100 ) {
    out << "CHARMM_PERIODIC " <<  func(x) << std::endl;
  }
  return Func::show_violations( out, x, verbose_level, threshold);

}


} //constraints
} //scoring
} //core
