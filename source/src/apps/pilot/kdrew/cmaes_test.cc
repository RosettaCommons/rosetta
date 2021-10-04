// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

// Utility headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

/*
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
*/


#include <stdio.h>
#include <stdlib.h> /* free() */
#include <stddef.h> /* NULL */
#include <cmaes/cmaes_interface.h>
#include <algorithm>

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

//using namespace core::conformation::symmetry;

// tracer - used to replace cout
static basic::Tracer TR( "CMAES_TEST" );


double fitfun(core::pose::Pose p, double const *x, int dim);

/* the objective (fitness) function to be minimized */
double fitfun(core::pose::Pose p, core::optimization::MinimizerMap min_map, double const *x, int N) { /* function "cigtab" */
	/* int i;
	double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
	for(i = 2; i < N; ++i)
	{
	printf("%f ", x[i]);
	sum += x[i]*x[i];
	}*/

	//kdrew: I might be mixing vectors of index 1 and 0
	core::optimization::Multivec dofs( min_map.nangles() );
	for ( int ii = 1; ii <= N; ++ii ) {
		dofs[ii] = x[ii-1];
		TR << dofs[ii] << ",";
	}
	TR << std::endl;

	//kdrew: clang build complains
	//TR << dofs << std::endl;

	min_map.copy_dofs_to_pose( p, dofs );
	core::scoring::ScoreFunctionOP score_fxn_ = scoring::get_score_function();
	TR << (*score_fxn_)(p) << std::endl;
	//printf("\n");
	return (*score_fxn_)(p);
}



int main(int argc, char *argv[]) {

	try {

		devel::init(argc, argv);
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,basic::options::option[basic::options::OptionKeys::in::file::s]()[1]);

		core::scoring::ScoreFunctionOP score_fxn_ = scoring::get_score_function();
		TR << (*score_fxn_)(pose) << std::endl;


		/*
		kinematics::MoveMapOP movemap = new kinematics::MoveMap();
		movemap->set_bb( true );
		movemap->set_chi( true );
		movemap->set_jump( false );

		assert( is_symmetric( pose.conformation() ) );
		*/

		core::optimization::MinimizerMap min_map;

		//kdrew: set up move_map
		core::kinematics::MoveMap move_map;
		move_map.set_bb(  true );
		move_map.set_jump( true );
		move_map.set_chi( true );

		min_map.setup( pose, move_map );
		score_fxn_->setup_for_minimizing( pose, min_map );
		core::optimization::Multivec dofs( min_map.nangles() );
		min_map.copy_dofs_from_pose( pose, dofs );
		TR << "nangles:" << min_map.nangles() << std::endl;
		TR << "dofs:"<< std::endl;
		//TR << dofs << std::endl;
		for ( int ii = 1; ii <= min_map.nangles(); ++ii ) {
			TR << dofs[ii] << ",";
		}
		TR << "end dofs"<< std::endl;


		cmaes_t evo; /* an CMA-ES type struct or "object" */
		double *arFunvals, *const*pop, *xfinal, *xbestever;
		int i;

		/* Initialize everything into the struct evo, 0 means default */
		/*arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "cmaes_initials.par");*/
		std::vector<double> xstart( min_map.nangles() );
		/*double rgsigma[1] = { 0.3 };*/
		std::vector<double> rgsigma( min_map.nangles() );
		std::fill (rgsigma.begin(), rgsigma.end(), 0.3);
		/*arFunvals = cmaes_init(&evo, 22, xstart, rgsigma, 0, 0, NULL);*/

		//kdrew: I might be mixing vectors of index 1 and 0
		for ( int ii = 1; ii <= min_map.nangles(); ++ii ) {
			xstart[ii-1] = dofs[ii];
		}

		for ( int ii = 0; ii < min_map.nangles(); ++ii ) {
			TR << xstart[ii] << ":";
			TR << rgsigma[ii] << ", ";

		}
		TR << std::endl;

		//kdrew: cmaes_init expects a double* but xstart and rgsigma are std::vectors, get pointer to first element
		arFunvals = cmaes_init(&evo, min_map.nangles(), &xstart[0], &rgsigma[0], 0, 0, NULL);

		TR << cmaes_SayHello(&evo) << std::endl;
		//cmaes_ReadSignals(&evo, "cmaes_signals.par");  /* write header and initial values */

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
			for ( i = 0; i < cmaes_Get(&evo, "lambda"); ++i ) {
				arFunvals[i] = fitfun(pose, min_map, pop[i], (int) cmaes_Get(&evo, "dim"));
			}

			/* update the search distribution used for cmaes_SamplePopulation() */
			cmaes_UpdateDistribution(&evo, arFunvals);

			/* read instructions for printing output or changing termination conditions */
			//cmaes_ReadSignals(&evo, "cmaes_signals.par");
			fflush(stdout); /* useful in MinGW */
		}
		TR << "Stop: " << cmaes_TestForTermination(&evo) << std::endl; /* print termination reason */

		//cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

		/* get best estimator for the optimum, xmean */
		xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
		for ( int ii = 1; ii <= min_map.nangles(); ++ii ) {
			dofs[ii] = xfinal[ii-1];
		}
		min_map.copy_dofs_to_pose( pose, dofs );
		std::stringstream xfinalpdbname;
		xfinalpdbname << "xfinal_pose.pdb";
		pose.dump_scored_pdb( xfinalpdbname.str(), *score_fxn_ );

		xbestever = cmaes_GetNew(&evo, "xbestever"); /* "xbestever" might be used as well */
		for ( int ii = 1; ii <= min_map.nangles(); ++ii ) {
			dofs[ii] = xbestever[ii-1];
		}
		min_map.copy_dofs_to_pose( pose, dofs );
		std::stringstream xbesteverpdbname;
		xbesteverpdbname << "xbestever_pose.pdb";
		pose.dump_scored_pdb( xbesteverpdbname.str(), *score_fxn_ );

		cmaes_exit(&evo); /* release memory */

		/* do something with final solution and finally release memory */
		free(xfinal);
		free(xbestever);


	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;

}
