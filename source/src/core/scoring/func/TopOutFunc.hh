// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/TopOutFunc.hh
/// @brief Implementation of phenix "top-out" function
///   Similar to Geman-McClure: harmonic near 'x0_', flat past 'limit_'
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_constraints_TopOutFunc_hh
#define INCLUDED_core_scoring_constraints_TopOutFunc_hh

#include <core/scoring/func/TopOutFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class TopOutFunc : public Func {
public:
	TopOutFunc( Real weight_in, Real x0_in, Real limit_in ) :
		x0_( x0_in ), weight_( weight_in), limit_( limit_in ) {}

	FuncOP
	clone() const { return new TopOutFunc( *this ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream & in );

	void show_definition( std::ostream &out ) const;

	Real x0() const { return x0_; }
	Real limit() const { return limit_; }
	Real weight() const { return weight_; }

	void x0( Real x ) { x0_ = x; }
	void limit( Real limit ) { limit_ = limit; }
	void weight( Real weight ) { weight_ = weight; }

	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	Real x0_;
	Real weight_;
	Real limit_;
};

} // constraints
} // scoring
} // core

#endif
