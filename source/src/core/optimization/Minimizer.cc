// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/Minimizer.cc
/// @brief  Minimizer class
/// @author Phil Bradley


// Unit headers
#include <core/optimization/LineMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/GA_Minimizer.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <utility/exit.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>


// C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>

#include <utility/vector1.hh>

#ifdef WIN32
#include <functional>
#endif


namespace core {
namespace optimization {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "core.optimization.Minimizer" );

// set the function and the options
Minimizer::Minimizer(
	Multifunc & func_in,
	MinimizerOptions const & options_in
) : func_( func_in ), options_( options_in ) {}

/////////////////////////////////////////////////////////////////////////////
/// See @ref minimization_overview "Minimization overview and concepts" for details.
Real
Minimizer::run(
	Multivec & phipsi_inout // starting position, and solution is returned here
) {
	// parse options
	std::string const type( options_.min_type() );

	Multivec phipsi( phipsi_inout ), dE_dphipsi( phipsi_inout );

	Real end_func;
	DFPMinConvergedFractional fractional_converge_test( options_.minimize_tolerance() );
	DFPMinConvergedAbsolute absolute_converge_test( options_.minimize_tolerance() );
	int const ITMAX( options_.max_iter() );

	if ( type == "linmin" ) {
		func_.dfunc( phipsi, dE_dphipsi );
		linmin( phipsi, dE_dphipsi, end_func, ITMAX );
	} else if ( type == "linmin_iterated" ) {
		linmin_iterated( phipsi, end_func, fractional_converge_test, ITMAX );
	} else if ( type == "linmin_iterated_atol" ) {
		linmin_iterated( phipsi, end_func, absolute_converge_test, ITMAX );
	} else if ( type == "dfpmin" ) {
		dfpmin( phipsi, end_func, fractional_converge_test, ITMAX );
	} else if ( type == "dfpmin_armijo" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, false, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		dfpmin_armijo( phipsi, end_func, fractional_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "dfpmin_armijo_nonmonotone" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, true, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		dfpmin_armijo( phipsi, end_func, fractional_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "dfpmin_atol" ) {
		dfpmin( phipsi, end_func, absolute_converge_test, ITMAX );
	} else if ( type == "dfpmin_armijo_atol" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, false, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		dfpmin_armijo( phipsi, end_func, absolute_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "dfpmin_armijo_nonmonotone_atol" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, true, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		dfpmin_armijo( phipsi, end_func, absolute_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "dfpmin_strong_wolfe" ) {
		LineMinimizationAlgorithmOP strong_wolfe_line_search( new StrongWolfeLineMinimization( func_, false, phipsi_inout.size() ) );
		dfpmin_armijo( phipsi, end_func, fractional_converge_test, strong_wolfe_line_search, ITMAX );
	} else if ( type == "dfpmin_strong_wolfe_atol" ) {
		LineMinimizationAlgorithmOP strong_wolfe_line_search( new StrongWolfeLineMinimization( func_, false, phipsi_inout.size() ) );
		dfpmin_armijo( phipsi, end_func, absolute_converge_test, strong_wolfe_line_search, ITMAX );
	} else if ( type == "lbfgs" ) {
		LineMinimizationAlgorithmOP brent_line_search( new BrentLineMinimization( func_, phipsi.size() ) );
		brent_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, fractional_converge_test, brent_line_search, ITMAX );
	} else if ( type == "lbfgs_armijo" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, false, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, fractional_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "lbfgs_armijo_rescored" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, true, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, fractional_converge_test, armijo_line_search, ITMAX, true );
	} else if ( type == "lbfgs_armijo_atol" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, false, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, absolute_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "lbfgs_armijo_nonmonotone" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, true, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, fractional_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "lbfgs_armijo_nonmonotone_atol" ) {
		LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( func_, true, phipsi_inout.size() ) );
		armijo_line_search->silent( options_.silent() );
		lbfgs( phipsi, end_func, absolute_converge_test, armijo_line_search, ITMAX );
	} else if ( type == "lbfgs_strong_wolfe" ) {
		LineMinimizationAlgorithmOP strong_wolfe_line_search( new StrongWolfeLineMinimization( func_, false, phipsi_inout.size() ) );
		lbfgs( phipsi, end_func, fractional_converge_test, strong_wolfe_line_search, ITMAX );
	} else if ( type == "GA" ) {
		GA_Minimizer gam(func_, options_);
		gam.run(phipsi, ITMAX);
	} else {
		utility_exit_with_message("unknown type of minimization '"+type+"'");
	}

	phipsi_inout = phipsi;
	return func_( phipsi );
}

////////////////////////////////////////////////////////////////////////
// Convergence test stuff
////////////////////////////////////////////////////////////////////////

bool
DFPMinConvergedFractional::operator()(
	Real Fnew,
	Real Fold
) {
	return ( 2.0f*std::abs( Fnew - Fold ) <=
		tolerance*( std::abs( Fnew ) + std::abs( Fold ) + eps )
	);
}

bool
DFPMinConvergedAbsolute::operator()(
	Real Fnew,
	Real Fold
) {
	return ( std::abs( Fnew - Fold ) <= tolerance );
}

/////////////////////////////////////////////////////////////////////////////
// Skeleton algorithm for multivariate optimization
/////////////////////////////////////////////////////////////////////////////
Real JJH_Minimizer::run(
	Multivec & current_position
) {
	int const problem_size( current_position.size() );
	int const max_iter( std::max( 200, problem_size/10 ));

	// reset storage variables for descent direction updates
	_get_direction.initialize();

	// Get the starting function value and gradient
	Multivec descent_direction( problem_size, 0.0);
	Real current_value( _func( current_position ) );
	_func.dfunc( current_position, descent_direction );
	// Convert to gradient (negative of the derivative) in-place
	std::transform( descent_direction.begin(), descent_direction.end(),
		descent_direction.begin(), std::negate<Real>() );

	// Iterate to convergence
	int iter( 0 );
	while ( iter++ < max_iter ) {
		Real previous_value = current_value;
		current_value = _line_min( current_position, descent_direction );

		if ( _converged(current_value, previous_value) ) return current_value;
		//   descent_direction = _get_direction();
	}

	TR.Warning << "WARNING: Minimization has exceeded " << max_iter << "iterations but has not converged!" << std::endl;
	return current_value;
}

// wrapper functions around minimization routines to allow for more flexibility
//void
//Minimizer::dfpmin(
//   Multivec & P,
//   Real & FRET,
//   ConvergenceTest & converge_test
//   ) const
//{
//  int const N( P.size() );
//  int const ITMAX = std::max( 200, N/10 );
//  Minimizer::dfpmin( P, FRET, converge_test, ITMAX );
//}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// now cut and paste from minimize.cc
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
void
Minimizer::dfpmin(
	Multivec & P,
	Real & FRET,
	ConvergenceTest & converge_test,
	int const ITMAX
) const {
	int const N( P.size() );
	// int const ITMAX = std::max( 200, N/10 );

	// Grab a line minimizer
	BrentLineMinimization test_brent( func_, N );
	LineMinimizationAlgorithm * line_min = &test_brent; // this variable is not neccessary

	// should get rid of these FArrays
	FArray2D< Real > HESSIN( N, N, 0.0 );
	Multivec XI( N );
	Multivec G( N );
	Multivec DG( N );

	// get function and its gradient
	Real FP;
	FP = func_(P);
	func_.dfunc(P,G);

	for ( int i = 1; i <= N; ++i ) {
		HESSIN(i,i) = 1.0;
		XI[i] = -G[i];
	}

	Multivec HDG( N );

	for ( int ITER = 1; ITER <= ITMAX; ++ITER ) {
		// Christophe added the following to allow premature end of minimization
		// I probably need to do the same with line_min
		if ( func_.abort_min(P) ) {
			TR.Warning << "WARNING: ABORTING MINIMIZATION TRIGGERED BY abort_min" << std::endl;
			return;
		}
		// End Christophe modifications

		// note that linmin modifes XI; afterward XI is the actual (vector)
		// step taken during linmin
		FRET = (*line_min)( P, XI );

		// check for convergence
		if ( converge_test( FRET, FP ) ) {
			//std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
			return;
		}

		FP = FRET;
		for ( int i = 1; i <= N; ++i ) {
			DG[i] = G[i];
		}

		// get function and its gradient
		FRET = func_(P);
		func_.dfunc(P,G);
		for ( int i = 1; i <= N; ++i ) {
			DG[i] = G[i]-DG[i];
		}
		for ( int i = 1; i <= N; ++i ) {
			HDG[i] = 0.0;
			for ( int j = 1; j <= N; ++j ) {
				HDG[i] += HESSIN(j,i)*DG[j];
			}
		}

		Real FAC, FAE, FAD;
		FAC = 0.;
		FAE = 0.;
		for ( int i = 1; i <= N; ++i ) {
			FAC += DG[i]*XI[i];
			FAE += DG[i]*HDG[i];
		}
		if ( FAC != 0.0 ) FAC = 1.0/FAC;
		if ( FAE != 0.0 ) {
			FAD = 1./FAE;
		} else {
			FAD = 0.;
		}
		for ( int i = 1; i <= N; ++i ) {
			DG[i] = FAC*XI[i] - FAD*HDG[i];
		}
		for ( int i = 1; i <= N; ++i ) {
			for ( int j = 1; j <= N; ++j ) {
				HESSIN(j,i) += FAC*XI[i]*XI[j] - FAD*HDG[i]*HDG[j] + FAE*DG[i]*DG[j];
			}
		}
		for ( int i = 1; i <= N; ++i ) {
			XI[i] = 0.;
			for ( int j = 1; j <= N; ++j ) {
				XI[i] -= HESSIN(j,i)*G[j];
			}
		}
	}

	if ( !options_.silent() ) {
		TR.Warning << "WARNING: DFPMIN MAX CYCLES " << ITMAX << " EXCEEDED, BUT FUNC NOT CONVERGED!" << std::endl;
	}

	//  std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
	return;

} // dfpmin
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// dfpmin Armijo
////////////////////////////////////////////////////////////////////////


void
Minimizer::dfpmin_armijo(
	Multivec & P,
	Real & FRET,
	ConvergenceTest & converge_test,
	LineMinimizationAlgorithmOP line_min,
	int const ITMAX
) const {
	int const N( P.size() );
	Real const EPS( 1.E-5 );

	FArray2D< Real > HESSIN( N, N, 0.0 );
	Multivec XI( N );
	Multivec G( N );
	Multivec DG( N );
	Multivec HDG( N );

	int const prior_func_memory_size( line_min->nonmonotone() ? 3 : 1 );
	Multivec prior_func_memory( prior_func_memory_size );

	if ( line_min->nonmonotone() ) line_min->_last_accepted_step = 0.005;

	// When inexact line search is used, HESSIN need not remain positive definite, so
	// additional safeguard must be added to ensure XI is a desc. direction (or inexact
	// line search would fail).  Two options for safeguards are implemented below:
	//   HOPT = 1  resets HESSIN to a multiple of identity when XI is not a desc. direction.
	//  HOPT = 2  leaves HESSIN unchanged if stepsize XMIN fails Wolfe's condition
	//         for ensuring new HESSIN to be positive definite.
	int const HOPT( 2 );

	// get function and its gradient
	int NF = 1;    // number of func evaluations
	Real prior_func_value = func_(P);
	func_.dfunc(P,G);

	// Start the prior function memory storage
	int func_memory_filled( 1 );
	prior_func_memory[ 1 ] = prior_func_value;

	for ( int i = 1; i <= N; ++i ) {
		HESSIN(i,i) = 1.0;
		XI[i] = -G[i];
	}

	Real FAC, FAE, FAD, FAF;

	for ( int ITER = 1; ITER <= ITMAX; ++ITER ) {
		line_min->_deriv_sum = 0.0;
		Real Gmax = 0.0;
		Real Gnorm = 0.0;

		for ( int i = 1; i <= N; ++i ) {
			line_min->_deriv_sum += XI[i]*G[i];
			Gnorm += G[i]*G[i];
			if ( std::abs( G[i] ) > Gmax ) {
				Gmax=std::abs( G[i] );
			}
		}

		Gnorm = std::sqrt(Gnorm);

		line_min->_func_to_beat = prior_func_memory[ 1 ];
		for ( int i = 2 ; i <= func_memory_filled ; ++i ) {
			if ( line_min->_func_to_beat < prior_func_memory[ i ] ) {
				line_min->_func_to_beat = prior_func_memory[ i ];
			}
		}

		// P is returned as new pt, and XI is returned as the change.
		FRET = (*line_min)( P, XI );

		// std::cout << "N= " << N << " ITER= " << ITER << " #F-eval= " << NF << " maxG= " << SS( Gmax ) << " Gnorm= " << SS( Gnorm ) << " step= " << SS( line_min->_last_accepted_step ) << " func= " << SS( FRET ) << std::endl;

		if ( converge_test( FRET, prior_func_value ) ) {
			//$$$   std::cout << "dfpmin called linmin " << linmin_count << " times" << std::endl;
			if ( Gmax<=options_.gmax_cutoff_for_convergence() ) {

				//std::cout << "N= " << N << " ITER= " << ITER << " #F-eval= " << NF << " maxG= " << SS( Gmax ) << " Gnorm= " << SS( Gnorm ) << " step= " << SS( line_min->_last_accepted_step ) << " func= " << SS( FRET ) << " time= " << SS( get_timer("dfpmin") ) << std::endl;

				//    std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
				return;
			} else {
				if ( std::abs(FRET-prior_func_value)<=EPS ) {
					Real XInorm = 0.0;
					for ( int i = 1; i <= N; ++i ) {
						XInorm += XI[i]*XI[i];
					}
					if ( line_min->_deriv_sum < -1e-3*Gnorm*XInorm ) {
						//      std::cout << "Failed line search while large _deriv_sum, quit! N= " << N << " ITER= " << ITER << " #F-eval= " << NF << " maxG= " << SS( Gmax ) << " Gnorm= " << SS( Gnorm ) << " step= " << SS( line_min->_last_accepted_step ) << " func= " << SS( FRET ) /*<< " time= " << SS( get_timer("dfpmin") )*/ << std::endl;

						//      std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
						return;
					}
					// Not convergence yet. Reinitialize HESSIN to a diagonal matrix & update direction XI.
					// This requires G to be correctly the gradient of the function.

					if ( !options_.silent() ) {
						TR.Warning << ":( reset HESSIN from failed line search" << std::endl;
					}

					line_min->_deriv_sum = 0.0;
					for ( int i = 1; i <= N; ++i ) {
						for ( int j = 1; j < i; ++j ) {
							HESSIN(j,i) = 0.0;
						}
						for ( int j = i+1; j <= N; ++j ) {
							HESSIN(j,i) = 0.0;
						}
						if ( HESSIN(i,i) < 0.01 ) HESSIN(i,i) = 0.01;
						if ( !options_.silent() ) {
							TR.Warning << "HESSIN for (i,i): " << i << ' ' << HESSIN(i,i) << std::endl;
							TR.Warning << "G for (i): " << i << ' ' << G[i] << std::endl;
						}
						// std::cerr << "HESSIN for (i,i): " << i << ' ' << HESSIN(i,i) << std::endl;
						// std::cerr << "G for (i): " << i << ' ' << G[i] << std::endl;
						XI[i] = -HESSIN(i,i)*G[i];
						line_min->_deriv_sum += XI[i]*G[i];
					}

					FRET = (*line_min)( P, XI );

					//     std::cout << "Failed line search again, quit! N= " << N << " ITER= " << ITER << " #F-eval= " << NF << " maxG= " << SS( Gmax ) << " Gnorm= " << SS( Gnorm ) << " step= " << SS( line_min->_last_accepted_step ) << " func= " << SS( FRET ) /*<< " time= " << SS( get_timer("dfpmin") )*/ << std::endl;

					if ( std::abs(FRET-prior_func_value)<=EPS ) {
						//      std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
						return;
					}
				}
			}
		}

		prior_func_value = FRET;

		// Update memory of function calls
		if ( func_memory_filled < prior_func_memory_size ) {
			func_memory_filled++;
		} else {
			for ( int i = 1 ; i < func_memory_filled ; ++ i ) {
				prior_func_memory[ i ] = prior_func_memory[ i + 1 ];
			}
		}
		prior_func_memory[ func_memory_filled ] = prior_func_value;

		for ( int i = 1; i <= N; ++i ) {
			DG[i] = G[i];
		}

		// Some line minimization algorithms require a curvature
		// check that involves the derivative before they accept a
		// move - in these cases we don't need to recalculate
		if ( line_min->provide_stored_derivatives() ) {
			line_min->fetch_stored_derivatives( G );
		} else {
			//FRET = func_(P);
			func_.dfunc(P,G);
		}

		NF++;

		line_min->_deriv_sum = 0.0;      //needed if HOPT = 2
		Real DRVNEW = 0.0;      //needed if HOPT = 2
		for ( int i = 1; i <= N; ++i ) {
			line_min->_deriv_sum += XI[i]*DG[i];   //needed if HOPT = 2
			DRVNEW += XI[i]*G[i];   //needed if HOPT = 2
			DG[i] = G[i]-DG[i];
		}

		//  if ( line_min->_last_accepted_step = 0.0 ) {
		//   std::cout << " line_min->_last_accepted_step = 0.0! " << std::endl; //diagnostic
		//  }

		if ( HOPT == 1 || DRVNEW > 0.95*line_min->_deriv_sum ) { //needed if HOPT = 2
			for ( int i = 1; i <= N; ++i ) {
				HDG[i] = 0.0;
				for ( int j = 1; j <= N; ++j ) {
					HDG[i] += HESSIN(j,i)*DG[j];
				}
			}
			FAC = 0.0;
			FAE = 0.0;
			FAF = 0.0;
			for ( int i = 1; i <= N; ++i ) {
				FAC += DG[i]*XI[i];
				FAE += DG[i]*HDG[i];
				FAF += DG[i]*DG[i];
			}
			FAF = FAC/FAF;
			FAC = 1.0/FAC;
			FAD = 1.0/FAE;
			for ( int i = 1; i <= N; ++i ) {
				DG[i] = FAC*XI[i] - FAD*HDG[i];
			}
			for ( int i = 1; i <= N; ++i ) {
				for ( int j = 1; j <= N; ++j ) {
					HESSIN(j,i) += FAC*XI[i]*XI[j] - FAD*HDG[i]*HDG[j] + FAE*DG[i]*DG[j];
				}
			}
		}         //needed if HOPT = 2

		for ( int i = 1; i <= N; ++i ) {
			XI[i] = 0.0;
			for ( int j = 1; j <= N; ++j ) {
				XI[i] -= HESSIN(j,i)*G[j];
			}
		}

		if ( HOPT == 1 ) {
			DRVNEW=0.0;
			for ( int i = 1; i <= N; ++i ) {
				DRVNEW += XI[i]*G[i];
			}
			// If direc. deriv >0, reset the Hessian inverse estimate
			if ( DRVNEW > -EPS ) {
				//    std::cout << "reset hessin; dirdg=" << SS( line_min->_deriv_sum ) << std::endl;
				if ( FAF<0.01 ) FAF=0.01;
				for ( int i = 1; i <= N; ++i ) {
					for ( int j = 1; j <= N; ++j ) {
						HESSIN(j,i) = 0;
					}
					HESSIN(i,i) = FAF;
					XI[i] = -FAF*G[i];
				}
			}
		} // HOPT == 1
	} // for ITER

	if ( !options_.silent() ) {
		TR.Warning << "WARNING: DFPMIN (Armijo) MAX CYCLES " << ITMAX << " EXCEEDED, BUT FUNC NOT CONVERGED!" << std::endl;
	}

	// std::cout << "Called line minimization " << line_min->_num_linemin_calls << std::endl;
	return;
}

////////////////////////////////////////////////////////////////////////
// *      Limited memory BFGS (L-BFGS).
// *
// * Copyright (c) 1990, Jorge Nocedal
// * Copyright (c) 2007-2010 Naoaki Okazaki
// * All rights reserved.
// *
// * Permission is hereby granted, free of charge, to any person obtaining a copy
// * of this software and associated documentation files (the "Software"), to deal
// * in the Software without restriction, including without limitation the rights
// * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// * copies of the Software, and to permit persons to whom the Software is
// * furnished to do so, subject to the following conditions:
// *
// * The above copyright notice and this permission notice shall be included in
// * all copies or substantial portions of the Software.
// *
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// * THE SOFTWARE.
////////////////////////////////////////////////////////////////////////

// This library is a C port of the FORTRAN implementation of Limited-memory
// Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
// The original FORTRAN source code is available at:
// http://www.ece.northwestern.edu/~nocedal/lbfgs.html
//
// The L-BFGS algorithm is described in:
//     - Jorge Nocedal.
//       Updating Quasi-Newton Matrices with Limited Storage.
//       <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
//     - Dong C. Liu and Jorge Nocedal.
//       On the limited memory BFGS method for large scale optimization.
//       <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

void
Minimizer::lbfgs(
	Multivec & X,
	Real & FRET,
	ConvergenceTest & converge_test,
	LineMinimizationAlgorithmOP line_min,
	int const ITMAX,
	bool w_rescore /*= false*/
) const {
	int const N( X.size() );
	static int M( basic::options::option[ basic::options::OptionKeys::optimization::lbfgs_M ]() );
	int const PAST( line_min->nonmonotone() ? 3 : 1 );
	Real const EPS( 1.E-5 );

	int K = 1; // number of func evaluations

	// Allocate working space.
	Multivec XP( N );
	Multivec Xtemp( N );
	Multivec G( N );
	Multivec GP( N );
	Multivec Gtemp( N );
	Multivec D( N );
	Multivec W( N );

	// Allocate & initialize limited memory storage
	int CURPOS = 1; // pointer to current location in lm
	utility::vector1< lbfgs_iteration_data > lm(M);
	for ( int i=1; i<=M; ++i ) {
		lm[i].alpha = 0;
		lm[i].ys = 0;
		lm[i].s.resize(N,0.0);
		lm[i].y.resize(N,0.0);
	}

	// Allocate space for storing previous values of the objective function
	Multivec pf(PAST);

	// Evaluate the function value and its gradient
	int func_memory_filled( 1 );
	Real prior_func_value = func_(X);
	pf[1] = FRET = prior_func_value;
	func_.dfunc(X,G);

	if ( w_rescore ) {  //fpd   reevaluate score after gradient computation
		pf[1] = FRET = prior_func_value = func_(X);
	}

	// Compute the direction
	// we assume the initial hessian matrix H_0 as the identity matrix.
	Real invdnorm = 0.0;
	for ( int i = 1; i <= N; ++i ) {
		D[i] = -G[i];
		//invdnorm += D[i]*D[i]; // Not used
	}
	//invdnorm = 1.0/sqrt( invdnorm ); // Not used

	if ( line_min->nonmonotone() ) line_min->_last_accepted_step = 0.005;

	// initial stepsize = 1/sqrt ( dot(D,D) )
	//line_min->_last_accepted_step = invdnorm;

	bool last_step_good = true;
	for ( int ITER = 1; ITER <= ITMAX; ++ITER ) {
		// Store the current position and gradient vectors
		if ( last_step_good ) {
			XP = X;
			GP = G;
		}

		// line min
		line_min->_deriv_sum = 0.0;
		Real Gmax = 0.0;
		//Real Gnorm = 0.0;
		for ( int i = 1; i <= N; ++i ) {
			line_min->_deriv_sum += D[i]*G[i];
			//Gnorm += G[i]*G[i]; //unused
			if ( std::abs( G[i] ) > Gmax ) {
				Gmax=std::abs( G[i] );
			}
		}
		//Gnorm = std::sqrt(Gnorm); //unused

		line_min->_func_to_beat = pf[ 1 ];
		for ( int i = 2 ; i <= func_memory_filled ; ++i ) {
			if ( line_min->_func_to_beat < pf[ i ] ) {
				line_min->_func_to_beat = pf[ i ];
			}
		}

		// check 1: if derivative is positive, flip signs of positive components
		if ( line_min->_deriv_sum > -EPS ) {
			//std::cerr << "reset 1" << std::endl;
			line_min->_deriv_sum = 0.0;
			for ( int i = 1; i <= N; ++i ) {
				if ( D[i]*G[i] >= 0 ) D[i] = -D[i];

				line_min->_deriv_sum += D[i]*G[i];
				//Gnorm += G[i]*G[i]; //unused
				if ( std::abs( G[i] ) > Gmax ) {
					Gmax=std::abs( G[i] );
				}
			}
			//Gnorm = std::sqrt(Gnorm); // unused
		}

		// check 2: if derivative still positive, reset Hessian
		if ( line_min->_deriv_sum > -EPS ) {
			//std::cerr << "reset 2" << std::endl;
			line_min->_deriv_sum = 0.0;
			for ( int i = 1; i <= N; ++i ) {
				D[i] = -G[i];
				line_min->_deriv_sum += D[i]*G[i];
			}

			if ( sqrt( -line_min->_deriv_sum ) > 1e-6 ) {
				//invdnorm = 1.0/sqrt( -line_min->_deriv_sum ); //unused
				//line_min->_last_accepted_step = invdnorm;  // keep prev stepsize
				func_memory_filled = 1;
				line_min->_func_to_beat = pf[1] = prior_func_value = FRET;
			}
		}


		// X is returned as new pt, and D is returned as the change
		FRET = (*line_min)( X, D );

		if ( converge_test( FRET, prior_func_value ) || line_min->_last_accepted_step == 0 ) {
			if ( Gmax<=options_.gmax_cutoff_for_convergence() ) {
				return;
			} else {
				if ( line_min->_last_accepted_step == 0 ) { // failed line search
					// Reset Hessian
					CURPOS = 1;
					K = 1;

					// reset line minimizer
					line_min->_deriv_sum = 0.0;
					for ( int i = 1; i <= N; ++i ) {
						D[i] = -G[i];
						line_min->_deriv_sum += D[i]*G[i];
					}

					if ( line_min->_deriv_sum <= -EPS ) {
						invdnorm = 1.0/sqrt( -line_min->_deriv_sum );

						// delete prior function memory
						line_min->_last_accepted_step = invdnorm; // reset initial step size (in case default of '1' is too large)
						func_memory_filled = 1;
						prior_func_value = FRET;
						pf[1] = prior_func_value;

						// line search in the direction of the gradient
						FRET = (*line_min)( X, D );
					} else {
						return;
					}

					// if the line minimzer fails again, abort
					if ( line_min->_last_accepted_step == 0 ) {
						if ( !options_.silent() ) {
							TR << "Line search failed even after resetting Hessian; aborting at iter#" << ITER << std::endl;
						}
						return;
					}
				}
			}
		}

		prior_func_value = FRET;

		// Update memory of function calls
		if ( func_memory_filled < PAST ) {
			func_memory_filled++;
		} else {
			for ( int i = 1 ; i < PAST ; ++ i ) {
				pf[ i ] = pf[ i + 1 ];
			}
		}
		pf[ func_memory_filled ] = prior_func_value;

		for ( int i = 1; i <= N; ++i ) {
			GP[i] = G[i];
		}

		// Some line minimization algorithms require a curvature
		// check that involves the derivative before they accept a
		// move - in these cases we don't need to recalculate
		if ( line_min->provide_stored_derivatives() ) {
			line_min->fetch_stored_derivatives( G );
		} else {
			func_.dfunc(X,G);
			if ( w_rescore ) {  //fpd   reevaluate score after gradient computation
				FRET = func_(X);
			}
		}

		line_min->_deriv_sum = 0.0;
		core::Real deriv_new = 0.0;
		for ( int i = 1; i <= N; ++i ) {
			line_min->_deriv_sum += D[i]*GP[i];
			deriv_new += D[i]*G[i];
		}

		// LBFGS updates
		//
		// Update vectors s and y:
		//   s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
		//   y_{k+1} = g_{k+1} - g_{k}.
		// Compute scalars ys and yy:
		//   ys = y^t \cdot s = 1 / \rho.
		//   yy = y^t \cdot y.
		// Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
		core::Real ys=0, yy=0;
		for ( int i = 1; i <= N; ++i ) {
			Xtemp[i] = X[i] - XP[i];
			Gtemp[i] = G[i] - GP[i];
			ys += Gtemp[i]*Xtemp[i];
			yy += Gtemp[i]*Gtemp[i];
		}

		if ( std::fabs( ys ) < 1e-6 ) {
			last_step_good = false;
		} else {
			last_step_good = true;
			for ( int i = 1; i <= N; ++i ) {
				lm[CURPOS].s[i] = Xtemp[i];
				lm[CURPOS].y[i] = Gtemp[i];
			}
			lm[CURPOS].ys = ys;

			// increment
			K++;
			CURPOS++;
			if ( CURPOS > M ) CURPOS = 1;
		}

		// Recursive formula to compute dir = -(H \cdot g).
		//   This is described in page 779 of:
		//   Jorge Nocedal.
		//   Updating Quasi-Newton Matrices with Limited Storage.
		//   Mathematics of Computation, Vol. 35, No. 151,
		//   pp. 773--782, 1980.
		int bound = std::min( M,K-1 );

		// Compute the negative of gradients
		for ( int i = 1; i <= N; ++i ) {
			D[i] = -G[i];
		}

		int j = CURPOS;
		for ( int pts=0; pts<bound; ++pts ) {
			j--;
			if ( j<=0 ) j=M; // wrap around
			//if (std::fabs(lm[j].ys) < 1e-6) continue;

			// \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}
			lm[j].alpha = 0;
			for ( int i = 1; i <= N; ++i ) {
				lm[j].alpha += lm[j].s[i] * D[i];
			}
			lm[j].alpha /= lm[j].ys;

			// q_{i} = q_{i+1} - \alpha_{i} y_{i}
			for ( int i = 1; i <= N; ++i ) {
				D[i] += -lm[j].alpha * lm[j].y[i];
			}
		}

		//for ( int i = 1; i <= N; ++i ) {
		// D[i] *= ys / yy;
		//}

		for ( int pts=0; pts<bound; ++pts ) {
			//if (std::fabs(lm[j].ys) < 1e-6) continue;

			// \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}
			core::Real beta=0.0;
			for ( int i = 1; i <= N; ++i ) {
				beta += lm[j].y[i] * D[i];
			}
			beta /= lm[j].ys;

			// \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}
			for ( int i = 1; i <= N; ++i ) {
				D[i] += (lm[j].alpha - beta) * lm[j].s[i];
			}

			j++;
			if ( j>M ) j=1; // wrap around
		}

	}

	if ( !options_.silent() ) {
		TR.Warning << "WARNING: LBFGS MAX CYCLES " << ITMAX << " EXCEEDED, BUT FUNC NOT CONVERGED!" << std::endl;
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////
void
Minimizer::linmin(
	Multivec & P,
	Multivec & XI,
	Real & FRET,
	int const
) const {
	//Try to use a line minimizer algorithm
	BrentLineMinimization test_brent( func_, P.size() );

	// See if this is good enough
	FRET = test_brent( P, XI );
}

////////////////////////////////////////////////////////////////////////
// linmin_iterated
////////////////////////////////////////////////////////////////////////

/// @brief Carry out many iterations of line searches, using only the gradient vector
/// (that is, don't approximate the Hessian or use it in any way).
/// @details Useful for testing minimization protocols.  In theory, this should converge
/// less efficiently than something that uses a Hessian approximation, but it will be entirely
/// free of any erratic behaviour that results from the approximation used.  Note that DFP/LBFGS
/// approximations are not guaranteed to work well with all functions.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
Minimizer::linmin_iterated(
	Multivec & P,
	Real & FRET,
	ConvergenceTest & converge_test,
	core::Size const ITMAX
) const {

	//Storage for the derivative:
	Multivec dE_dphipsi(P.size(), 0.0);

	//The line minimizer:
	BrentLineMinimization test_brent( func_, P.size() );
	test_brent.set_deriv_cutoff( options_.linmin_deriv_cutoff() );

	core::Size ITER(1);

	for ( ITER = 1; ITER <= ITMAX; ++ITER ) {
		core::Real const prior_func_value = func_(P);
		func_.dfunc( P, dE_dphipsi );

		//Perform line minimization and return final value:
		FRET = test_brent( P, dE_dphipsi );

		// Check for convergence
		if ( converge_test( FRET, prior_func_value ) ) {
			break;
		}
	}

	if ( ITER == ITMAX ) {
		TR.Warning << "Warning: maximum number of minimization iterations (" << ITMAX << ") was reached.  Function likely not converged." << std::endl;
		TR.Warning.flush();
	}

	return;
} //linmin_iterated

////////////////////////////////////////////////////////////////////////


} // namespace optimization
} // namespace core
