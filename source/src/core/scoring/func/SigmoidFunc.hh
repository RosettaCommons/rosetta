// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SigmoidFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author

#ifndef INCLUDED_core_scoring_func_SigmoidFunc_hh
#define INCLUDED_core_scoring_func_SigmoidFunc_hh

#include <core/scoring/func/SigmoidFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

class SigmoidFunc : public Func {
public:
	SigmoidFunc( Real const x0_in, Real const slope_in ): x0_( x0_in ), slope_( slope_in ){}

	FuncOP
	clone() const { return FuncOP( new SigmoidFunc( *this ) ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream & in );

	Size show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const;
	void show_definition( std::ostream &out ) const;

	Real x0() const {
		return x0_;
	}

	Real slope() const {
		return slope_;
	}

	void x0( Real x ) {
		x0_ = x;
	}

	void slope( Real slope ) {
		slope_ = slope;
	}

private:
	Real x0_;
	Real slope_;
};

} // constraints
} // scoring
} // core

#endif
