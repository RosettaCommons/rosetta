// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/SigmoidalFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author RobertVernon

#ifndef INCLUDED_core_scoring_func_CircularSigmoidalFunc_hh
#define INCLUDED_core_scoring_func_CircularSigmoidalFunc_hh

#include <core/scoring/func/CircularSigmoidalFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

/// @brief Function that operates in radians, for applications like DihedralConstraint.
/// Prevents discontinuities at 0/360 or -180/180 degrees for dihedral constraints.
class CircularSigmoidalFunc : public Func {
public:
	CircularSigmoidalFunc(
		Real const center_radians, Real const width_radians, Real const slope_radians
	) : xC_( center_radians ), m_( slope_radians ),
		o1_( -width_radians/2 ), o2_( width_radians/2 ), offset_( 1.0 ) {}

	CircularSigmoidalFunc(
		Real const center_radians, Real const width_radians, Real const slope_radians,
		Real const offset         ) : xC_( center_radians ), m_( slope_radians ),
		o1_( -width_radians/2 ), o2_( width_radians/2 ), offset_( offset ) {}

	FuncOP clone() const { return FuncOP( new CircularSigmoidalFunc( *this ) ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	virtual void read_data( std::istream & in );
	virtual void show_definition( std::ostream & out ) const;

private:
	Real xC_;
	Real m_;
	Real o1_;
	Real o2_;
	Real offset_;
};

} // constraints
} // scoring
} // core

#endif
