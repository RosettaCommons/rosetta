// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/IdentityFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_constraints_IdentityFunc_hh
#define INCLUDED_core_scoring_constraints_IdentityFunc_hh

#include <core/scoring/func/IdentityFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class IdentityFunc : public Func {
public:
	IdentityFunc() {}

	FuncOP
	clone() const {
		return new IdentityFunc( *this );
	}

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream & in );

	void show_definition( std::ostream &out ) const;

	Size
	show_violations(
		std::ostream & out, Real x, Size verbose_level, core::Real threshold = 1
	) const;

private:
	Real x0_;
	Real sd_;
};

} // constraints
} // scoring
} // core

#endif
