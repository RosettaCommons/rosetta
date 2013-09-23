// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/GaussianChainSingleFunc.hh
/// @brief Definition for functions used in definition of loop closure energies.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_constraints_GaussianChainSingleFunc_HH
#define INCLUDED_core_scoring_constraints_GaussianChainSingleFunc_HH

#include <core/scoring/constraints/GaussianChainSingleFunc.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// See GaussianChainFunc.cc for more information, including link to mathematical derivation.

// C++ Headers
namespace core {
namespace scoring {
namespace constraints {

class GaussianChainSingleFunc : public Func {
public:

	GaussianChainSingleFunc( Real const gaussian_variance_ );

	FuncOP
	clone() const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

private:

	void initialize_parameters();
	void recompute_parameters();

private:

	Real gaussian_variance_;

	// following have nice default values
	Real kB_T_;
	Real loop_fixed_cost_;

	// derived from above.
	Real loop_fixed_cost_total_;

};



} // constraints
} // scoring
} // core

#endif
