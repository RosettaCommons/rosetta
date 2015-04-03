// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/GaussianChainFunc.hh
/// @brief Definition for functions used in definition of loop closure energies.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_func_GaussianChainFunc_HH
#define INCLUDED_core_scoring_func_GaussianChainFunc_HH

#include <core/scoring/func/GaussianChainFunc.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// See GaussianChainFunc.cc for more information, including link to mathematical derivation.

// C++ Headers
namespace core {
namespace scoring {
namespace func {

class GaussianChainFunc : public Func {
public:

	GaussianChainFunc(
			 Real const gaussian_variance,
			 Real const loop_fixed_cost,
			 utility::vector1< Real > const & other_distances );

	GaussianChainFunc( Real const gaussian_variance,
										 Real const loop_fixed_cost );

	FuncOP
	clone() const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

	void set_force_combined_gaussian_approximation( bool const & setting ){ force_combined_gaussian_approximation_ = setting; initialize_func(); }

	bool force_combined_gaussian_approximation() const{ return force_combined_gaussian_approximation_; }

private:

	void initialize_parameters();
	void initialize_func();

private:

	// needed to define function.
	Real gaussian_variance_;
	Real loop_fixed_cost_;
	utility::vector1< Real > other_distances_;

	// following have nice default values.
	Real kB_T_;
	bool force_combined_gaussian_approximation_;

	FuncOP func_;
};


} // constraints
} // scoring
} // core

#endif
