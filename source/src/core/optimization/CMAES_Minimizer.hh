// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.



/// @file   core/optimization/CMAES_Minimizer.hh
/// @brief  Minimizer based on Covariance Matrix Adaptation Evolution Strategy (CMAES)
/// @author Kevin Drew

#ifndef INCLUDED_core_optimization_CMAES_Minimizer_hh
#define INCLUDED_core_optimization_CMAES_Minimizer_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/optimization/Multifunc.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {


class CMAES_Minimizer
{
public:
	CMAES_Minimizer(Multifunc & func_in, MinimizerOptions const & options):
		func_( func_in ),
		minimize_tolerance_( options.minimize_tolerance() ),
		rgsigma_(options.cmaes_rgsigma()),
		max_iter_(options.max_iter())

	{}
	/// @brief run full blackbox CMAES minimization and return final score.
	//  CMAES minimizer algorithm will sample a population around the starting conformation and after each
	//  iteration will update the distribution from which the samples are drawn. The update is weighted
	//  towards lower energy samples so as to converge to a minimum.
	//  The outer while loop tests for convergence and the inner for loop samples the distribution.
	Real run( Multivec & phipsi_inout); // starting position, and solution is returned here

	void rgsigma( Real rgsigma_in );

private:
	Multifunc & func_;

	Real minimize_tolerance_;
	Real rgsigma_;
	Real max_iter_;

};


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_CMAES_Minimizer_hh
