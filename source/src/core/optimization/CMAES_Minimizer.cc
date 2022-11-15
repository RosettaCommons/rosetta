// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/CMAES_Minimizer.cc
/// @brief  Minimizer based on Covariance Matrix Adaptation Evolution Strategy (CMAES)
/// @author Kevin Drew


#include <core/optimization/CMAES_Minimizer.hh>

#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <algorithm>

#include <core/optimization/Multifunc.hh>
#include <utility/vector1.hh>

#include <cmaes/cmaes_interface.h>


namespace core {
namespace optimization {

static basic::Tracer TR( "core.optimization.CMAES_Minimizer" );


using core::Size;


void CMAES_Minimizer::rgsigma( Real rgsigma_in )
{
	rgsigma_=rgsigma_in;
}

void CMAES_Minimizer::lambda( int lambda_in )
{
	lambda_=lambda_in;
}


// starting position, and solution is returned here
Real CMAES_Minimizer::run( Multivec & v )
{
	cmaes_t evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xfinal;


	/* Initialize everything into the struct evo, 0 means default */
	std::vector<double> xstart( v.size() );
	std::vector<double> rgsigma( v.size() );
	std::fill (rgsigma.begin(), rgsigma.end(), rgsigma_);

	//kdrew: Carefull mixing vectors of index 1 and 0
	for ( Size ii = 1; ii <= v.size(); ++ii ) {
		xstart[ii-1] = v[ii];
	}

	/*for ( int ii = 0; ii < min_map.nangles(); ++ii )
	{
	std::cout << xstart[ii] << ":";
	std::cout << rgsigma[ii] << ", ";

	}
	std::cout << std::endl;*/

	//kdrew: interface expects a parameter filename but we aren't use it, clang complains about NULL pointers so just give it an empty filename for now
	const char *input_parameter_filename = "";

	//kdrew: cmaes_init expects a double* but xstart and rgsigma are std::vectors, get pointer to first element
	//kdrew: cmaes_init looks for cmaes_initials.par by default and if it doesn't find anything, it adds an error line to errcmaes.err, ignoring for now
	//kdrew: lambda default is 0 but will be converted to 4+(int)(3*log((double)N)), where N is size of dims in the c-cmaes library code
	arFunvals = cmaes_init(&evo, v.size(), &xstart[0], &rgsigma[0], numeric::random::rg().random_range(0,9999999), lambda_, input_parameter_filename);

	TR << "CMAES_HELLO: " << cmaes_SayHello(&evo) << std::endl;
	//cmaes_ReadSignals(&evo, "cmaes_signals.par");  /* write header and initial values */

	//kdrew: set termination conditions
	evo.sp.stopMaxIter = max_iter_; //kdrew: default is 1e299
	evo.sp.stopTolFun = minimize_tolerance_;

	TR << "rgsigma: " << rgsigma_ << std::endl;
	TR << "lambda: " << lambda_ << std::endl;
	//TR << "stopMaxIter: " << stopMaxIter_ << std::endl;
	TR << "stopMaxIter: " << max_iter_ << std::endl;
	TR << "stopTolFun: " << minimize_tolerance_ << std::endl;


	/* Iterate until stop criterion holds */
	while ( !cmaes_TestForTermination(&evo) )
			{
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

		/* Here we may resample each solution point pop[i] until it
		becomes feasible. function is_feasible(...) needs to be
		user-defined.
		Assumptions: the feasible domain is convex, the optimum is
		not on (or very close to) the domain boundary, initialX is
		feasible and initialStandardDeviations are sufficiently small
		to prevent quasi-infinite looping. */
		/* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
		while (!is_feasible(pop[i]))
		cmaes_ReSampleSingle(&evo, i);
		*/

		/* evaluate the new search points using fitfun */
		//for ( int i = 0; i < cmaes_Get(&evo, "lambda"); ++i ) {
		for ( int i = 0, imax = cmaes_Get(&evo, "lambda"); i < imax; ++i ) {

			core::optimization::Multivec vtrial( (int) cmaes_Get(&evo, "dim") );
			for ( int ii = 1; ii <= (int) cmaes_Get(&evo, "dim"); ++ii ) {
				vtrial[ii] = pop[i][ii-1];
			}

			arFunvals[i] = func_( vtrial );
		}

		/* update the search distribution used for cmaes_SamplePopulation() */
		cmaes_UpdateDistribution(&evo, arFunvals);

		/* read instructions for printing output or changing termination conditions */
		//cmaes_ReadSignals(&evo, "cmaes_signals.par");
		fflush(stdout); /* useful in MinGW */
	}
	//printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
	TR << "Stop: " << cmaes_TestForTermination(&evo) << std::endl; /* print termination reason */

	//cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

	/* get best estimator for the optimum, xmean */
	xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
	for ( Size ii = 1; ii <= v.size(); ++ii ) {
		v[ii] = xfinal[ii-1];
	}

	return func_(v);
}


} // namespace optimization
} // namespace core
