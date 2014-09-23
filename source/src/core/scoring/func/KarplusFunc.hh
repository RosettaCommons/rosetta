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
/// @author Nikolas Sgourakis

#ifndef INCLUDED_core_scoring_func_KarplusFunc_hh
#define INCLUDED_core_scoring_func_KarplusFunc_hh

#include <core/scoring/func/KarplusFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

/// @brief Function that evaluates a J-coupling from dihedral angles in radians, for applications like DihedralConstraint.
/// Prevents discontinuities at 0/360 or -180/180 degrees for dihedral constraints.
class KarplusFunc : public Func {
public:

	KarplusFunc(
	  Real const A_Hertz , Real const B_Hertz, Real const C_Hertz, Real const Dphi_radians, Real const x0_Hertz, Real const sd_Hertz, Real const offset=0.0
											 ):A_(A_Hertz), B_(B_Hertz), C_(C_Hertz), Dphi_(Dphi_radians), x0_( x0_Hertz ), sd_( sd_Hertz ), offset_( offset ) {}

	FuncOP clone() const { return FuncOP( new KarplusFunc( *this ) ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	virtual void read_data( std::istream & in );
	virtual void show_definition( std::ostream & out ) const;

private:
  Real A_;
  Real B_;
  Real C_;
  Real Dphi_;
	Real x0_;
	Real sd_;
	Real offset_;
};

} // constraints
} // scoring
} // core

#endif
