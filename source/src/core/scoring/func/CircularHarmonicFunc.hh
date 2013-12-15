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

#ifndef INCLUDED_core_scoring_func_CircularHarmonicFunc_hh
#define INCLUDED_core_scoring_func_CircularHarmonicFunc_hh

#include <core/scoring/func/CircularHarmonicFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

/// @brief Function that operates in radians, for applications like DihedralConstraint.
/// Prevents discontinuities at 0/360 or -180/180 degrees for dihedral constraints.
class CircularHarmonicFunc : public Func {
public:
	CircularHarmonicFunc(
		Real const x0_radians, Real const sd_radians
											 ): x0_( x0_radians ), sd_( sd_radians ), offset_( 0.0 ) {}

	CircularHarmonicFunc(
											 Real const x0_radians, Real const sd_radians, Real const offset
											 ): x0_( x0_radians ), sd_( sd_radians ), offset_( offset ) {}

	FuncOP clone() const { return new CircularHarmonicFunc( *this ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	virtual void read_data( std::istream & in );
	virtual void show_definition( std::ostream & out ) const;

	Real x0() const {
		return x0_;
	}

	Real sd() const {
		return sd_;
	}

private:
	Real x0_;
	Real sd_;
	Real offset_;
};

} // constraints
} // scoring
} // core

#endif
