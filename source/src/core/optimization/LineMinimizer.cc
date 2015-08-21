// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/LineMinimizer.cc
/// @brief  line minimizer classes
/// @author Phil Bradley
/// @author Jim Havranek


// Unit headers
#include <core/optimization/LineMinimizer.hh>

#include <ObjexxFCL/Fmath.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <cmath>
//#include <cstdlib>
// #include <cstdio>
#include <iostream>
#include <algorithm>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL {
} using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "core.optimization.LineMinimizer" );

namespace core {
namespace optimization {

/// @details Auto-generated virtual destructor
LineMinimizationAlgorithm::~LineMinimizationAlgorithm() {}

void
func_1d::dump( Real displacement ) {
	for ( uint i =  1 ; i <= _starting_point.size() ; ++i ) {
		_eval_point[i] = _starting_point[i] + ( displacement * _search_direction[i] );
	}
	return _func.dump( _starting_point, _eval_point );
}

static thread_local basic::Tracer TR( "core.optimization.LineMinimizer" );


// Functor'ed up version of accurate line minimization
/////////////////////////////////////////////////////////////////////////////
Real
BrentLineMinimization::operator()(
	Multivec & current_position,
	Multivec & search_direction
)
{
	int const problem_size( current_position.size() );

	_num_linemin_calls++;

	Real AX, XX, BX, FA, FX, FB;

	// check magnitude of derivative
	Real derivmax = 0.0;
	for ( int i = 1; i <= problem_size; ++i ) {
		if ( std::abs(search_direction[i]) > std::abs(derivmax) ) {
			derivmax = search_direction[i];
		}
	}

	if ( std::abs(derivmax) <= .0001 ) {
		Real final_value =_func(current_position);
		return final_value; // deriv = 0, return value
	}

	// Construct the one-dimensional projection of the function
	func_1d this_line_func( current_position, search_direction, _func );

	// initial step sizes from our options
	AX = _ax;
	XX = _xx; // initial range for bracketing
	BX = _bx; // now set in namespace in  minimize_set_func

	MNBRAK(AX,XX,BX,FA,FX,FB,this_line_func);

	Real final_value = BRENT(AX,XX,BX,FA,FX,FB,_tolerance,this_line_func);

	for ( int j = 1; j <= problem_size; ++j ) {
		search_direction[j] *= _last_accepted_step;
		current_position[j] += search_direction[j];
	}

	//  TR.Info << "Linemin used " << this_line_func.get_eval_count() <<
	//    " function calls" << std::endl;

	return final_value;
}

/////////////////////////////////////////////////////////////////////////////
void
BrentLineMinimization::MNBRAK(
	Real & AX,
	Real & BX,
	Real & CX,
	Real & FA,
	Real & FB,
	Real & FC,
	func_1d & func_eval
) const
{
	Real DUM, R, Q, U, ULIM, FU;
	Real const GOLD = { 1.618034 };
	Real const GLIMIT = { 100.0 };
	Real const TINY = { 1.E-20 };

	FA = func_eval(AX);
	FB = func_eval(BX);
	if ( FB > FA ) {
		DUM = AX;
		AX = BX;
		BX = DUM;
		DUM = FB;
		FB = FA;
		FA = DUM;
	}
	CX = BX+GOLD*(BX-AX);
	FC = func_eval(CX);

	while ( FB >= FC ) {
		R = (BX-AX)*(FB-FC);
		Q = (BX-CX)*(FB-FA);
		U = BX-((BX-CX)*Q-(BX-AX)*R)/(2.*sign(std::max(std::abs(Q-R),TINY),Q-R));
		ULIM = BX+GLIMIT*(CX-BX);
		if ( (BX-U)*(U-CX) > 0.0 ) {
			FU = func_eval(U);
			if ( FU < FC ) {
				AX = BX;
				FA = FB;
				BX = U;
				FB = FU;
				continue;
			} else if ( FU > FB ) {
				CX = U;
				FC = FU;
				continue;
			}
			U = CX+GOLD*(CX-BX);
			FU = func_eval(U);

		} else if ( (CX-U)*(U-ULIM) > 0. ) {
			FU = func_eval(U);
			if ( FU < FC ) {
				BX = CX;
				CX = U;
				U = CX+GOLD*(CX-BX);
				FB = FC;
				FC = FU;
				FU = func_eval(U);
			}
		} else if ( (U-ULIM)*(ULIM-CX) >= 0. ) {
			U = ULIM;
			FU = func_eval(U);
		} else {
			U = CX+GOLD*(CX-BX);
			FU = func_eval(U);
		}
		AX = BX;
		BX = CX;
		CX = U;
		FA = FB;
		FB = FC;
		FC = FU;
	}
} // mnbrak

/////////////////////////////////////////////////////////////////////////////
//
Real
BrentLineMinimization::BRENT(
	Real const AX,
	Real const BX,
	Real const CX,
	Real & FA,
	Real & FB,
	Real const FC,
	Real const TOL,
	func_1d & func_eval
)
{
	Real const brent_abs_tolerance( _abs_tolerance );

	Real BRENT; // Return value

	Real A, B, E,TOL1,TOL2;
	Real V,W,X,FX,FV,FW,XM,R,Q,P,ETEMP,D,U,FU;

	int const ITMAX = { 100 };
	Real const CGOLD = { 0.3819660 };
	Real const ZEPS = { 1.0E-10 };
	//$$$      int func_counter;       // diagnostic only

	A = std::min(AX,CX);
	B = std::max(AX,CX);
	V = BX;
	W = V;
	X = V;
	E = 0.0;
	//     these two initializations added to remove warnings
	D = 0.0; // D will be set first time through loop, when if (E>tol) is false
	BRENT = -999.9;

	//********************************************
	FX = FB;
	if ( A == AX ) {
		FB = FC;
	} else {
		FB = FA;
		FA = FC;
	}

	// pb 3/20/07: this was already commented out:
	//      FX = F(gfrag,X);
	//********************************************
	//$$$      func_counter = 1;

	FV = FX;
	FW = FX;
	for ( int iter = 1; iter <= ITMAX; ++iter ) {
		XM = 0.5*(A+B);
		TOL1 = TOL*std::abs(X)+ZEPS;
		TOL2 = 2.*TOL1;


		//********************************************
		//     here, we exit BRENT if the function varies by less than a
		//     tolerance over the interval
		if ( ( std::abs(FB-FX) <= brent_abs_tolerance ) &&
				( std::abs(FA-FX) <= brent_abs_tolerance ) ) break;

		//********************************************
		//     here, we exit BRENT if we have narrowed the search down to a
		//     small change in X
		if ( std::abs(X-XM) <= (TOL2-.5*(B-A)) ) break;

		bool skipnow = false; //To replace a goto statement that was being used to skip a few lines of code.
		if ( std::abs(E) > TOL1 ) {
			R = (X-W)*(FX-FV);
			Q = (X-V)*(FX-FW);
			P = (X-V)*Q-(X-W)*R;
			Q = 2.*(Q-R);
			if ( Q > 0.0 ) P = -P;
			Q = std::abs(Q);
			ETEMP = E;
			E = D;
			if ( ! ( std::abs(P) >= std::abs(.5*Q*ETEMP) ||
					P <= Q*(A-X) || P >= Q*(B-X) ) ) {
				D = P/Q;
				U = X+D;
				if ( U-A < TOL2 || B-U < TOL2 ) D = sign(TOL1,XM-X);
				skipnow = true;
			}
		}
		if ( !skipnow ) {
			if ( X >= XM ) {
				E = A-X;
			} else {
				E = B-X;
			}
			D = CGOLD*E;
		}
		if ( std::abs(D) >= TOL1 ) {
			U = X+D;
		} else {
			U = X+sign(TOL1,D);
		}

		// call Minimizer object's 1D function using stored data
		FU = func_eval(U);
		//$$$        ++func_counter;

		if ( FU <= FX ) {
			if ( U >= X ) {
				A = X;
				FA = FX;
			} else {
				B = X;
				FB = FX;
			}
			V = W;
			FV = FW;
			W = X;
			FW = FX;
			X = U;
			FX = FU;
		} else {
			if ( U < X ) {
				A = U;
				FA = FU;
			} else {
				B = U;
				FB = FU;
			}
			if ( FU <= FW || W == X ) {
				V = W;
				FV = FW;
				W = U;
				FW = FU;
			} else if ( FU <= FV || V == X || V == W ) {
				V = U;
				FV = FU;
			}
		}

		if ( iter >= ITMAX ) {
			TR.Error << "BRENT exceed maximum iterations. " << iter << std::endl;
			std::exit( EXIT_FAILURE );
		}
	} // iter=1,ITMAX


	_last_accepted_step = X;
	BRENT = FX;

	return BRENT;
}

/////////////////////////////////////////////////////////////////////////////
//
Real
ArmijoLineMinimization::operator()(
	Multivec & current_position,
	Multivec & search_direction
)
{

	Real const FACTOR( 0.5 );
	int const problem_size( current_position.size() );
	// Real max_step_limit( _nonmonotone ? 2.0 : 1.0 );
	Real max_step_limit( 1.0 );

	_num_linemin_calls++;

	// Construct the one-dimensional projection of the function
	func_1d this_line_func( current_position, search_direction, _func );

	// Early termination for derivatives (search_direction magnitudes) near zero
	// Please note that the search_direction vector is not normalized
	Real derivmax = 0.0;
	for ( int i = 1 ; i <= problem_size; ++i ) {
		if ( std::abs( search_direction[ i ] ) >
				std::abs( derivmax ) ) derivmax = search_direction[ i ];
	}

	// if ( runlevel >= gush) std::cout << "derivmax," << SS( derivmax ) << std::endl;
	if ( std::abs(derivmax) < .0001 ) {
		Real final_value = this_line_func( 0.0 ); // deriv = 0, return value
		return final_value;
	}

	//initial trial stepsize
	Real init_step( _last_accepted_step / FACTOR );
	if ( init_step > max_step_limit ) init_step = max_step_limit;

	Real final_value = Armijo( init_step, this_line_func );

	for ( int j = 1; j <= problem_size; ++j ) {
		search_direction[j] *= _last_accepted_step;
		current_position[j] += search_direction[j];
	}

	// std::cout << "Linemin used " << this_line_func.get_eval_count() <<
	//  " function calls and returns " << final_value << " on step of " << _last_accepted_step << std::endl;

	return final_value;
}

/////////////////////////////////////////////////////////////////////////////
//
Real
ArmijoLineMinimization::Armijo(
	Real init_step,
	func_1d & func_eval
)
{
	// INPUT PARAMETERS: XX,func,gfrag,NF
	// OUTPUT PARAMETERS: NF,FRET
	//
	// Given a function FUNC, and initial stepsize XX
	// such that 0 < XX and FUNC has negative derivative _deriv_sum at 0,
	// this routine returns a stepsize XX at least as good as that
	// given by Armijo rule.
	// Reference:  D.P. Bertsekas, Nonlinear Programming, 2nd ed, 1999, page 29.
	Real const FACTOR( 0.5 );
	Real const SIGMA( 0.1 );
	Real const SIGMA2( 0.8 );
	static Real const MINSTEP( basic::options::option[ basic::options::OptionKeys::optimization::armijo_min_stepsize ]() );  // default 1e-8

	//std::cout << "func_to_beat is " << _func_to_beat << std::endl;

	Real func_value = func_eval( init_step );
	_num_calls++;

	_last_accepted_step = init_step;

	if ( func_value < _func_to_beat+init_step*SIGMA2*_deriv_sum ) {
		Real test_step = init_step/FACTOR;
		Real test_func_value = func_eval( test_step );
		_num_calls++;
		if ( test_func_value < func_value ) {
			_last_accepted_step = test_step;
			return test_func_value;
		}
		return func_value;
	}

	Real far_step = init_step;
	while ( func_value > _func_to_beat + init_step*SIGMA*_deriv_sum ) {
		// Abort if function value is unlikely to improve.
		if ( ( init_step <= 1e-5 * far_step ) ||
				(init_step < MINSTEP && func_value >= _func_to_beat) ) {
			Real test_step = ( func_value - _func_to_beat ) / init_step;
			if ( !_silent ) {
				TR.Error << TR.Red << "Inaccurate G! step= " << ( init_step ) << " Deriv= " <<
					( _deriv_sum ) << " Finite Diff= " << ( test_step ) << TR.Reset << std::endl;
			}
			func_eval.dump( init_step );
			_last_accepted_step = 0.0;
			return _func_to_beat;
		}

		init_step *= FACTOR*FACTOR;  // faster step decrease
		// init_step *= FACTOR;
		//std::cout << "func_eval( " << init_step << ")" << std::endl;
		func_value = func_eval( init_step );
		_num_calls++;
	}

	_last_accepted_step = init_step;

	if ( init_step < 0.0 ) {
		TR << "Forced to do parabolic fit!" << std::endl;
		// Parabola interpolate between 0 and init_step for refinement
		Real test_step = -_deriv_sum*init_step*init_step/
			(2*(func_value - _func_to_beat - init_step * _deriv_sum));
		if ( test_step > 1e-3*far_step && test_step < far_step ) {
			Real test_func_value = func_eval( test_step );
			_num_calls++;
			if ( test_func_value < func_value ) {
				_last_accepted_step = test_step;
				func_value = test_func_value;
			}
		}
	}

	return func_value;
}

/////////////////////////////////////////////////////////////////////////////
//
Real
StrongWolfeLineMinimization::operator()(
	Multivec & current_position,
	Multivec & search_direction
)
{

	int const problem_size( current_position.size() );
	// Real max_step_limit( _nonmonotone ? 2.0 : 1.0 );

	_num_linemin_calls++;

	// Construct the one-dimensional projection of the function
	func_1d this_line_func( current_position, search_direction, _func );

	// Early termination for derivatives (search_direction magnitudes) near zero
	// Please note that the search_direction vector is not normalized
	Real derivmax = 0.0;
	for ( int i = 1 ; i <= problem_size; ++i ) {
		if ( std::abs( search_direction[ i ] ) >
				std::abs( derivmax ) ) derivmax = search_direction[ i ];
	}

	// if ( runlevel >= gush) std::cout << "derivmax," << SS( derivmax ) << std::endl;
	if ( std::abs(derivmax) < .0001 ) {
		Real final_value = this_line_func( 0.0 ); // deriv = 0, return value
		//  TR << "Exiting line minimization due to super small derivative" << std::endl;
		return final_value;
	}

	//initial trial stepsize
	Real init_step( std::min( 2.0*_last_accepted_step, 1.0 ) );

	Real final_value = StrongWolfe( init_step, this_line_func );

	for ( int j = 1; j <= problem_size; ++j ) {
		search_direction[j] *= _last_accepted_step;
		current_position[j] += search_direction[j];
	}

	// TR << "Linemin used " << this_line_func.get_eval_count() <<
	//  " function calls and " << this_line_func.get_deriv_count() << " deriv calls and  returns " << final_value <<
	//  " on step of " << _last_accepted_step << std::endl;

	return final_value;
}

/////////////////////////////////////////////////////////////////////////////
//
Real
StrongWolfeLineMinimization::StrongWolfe(
	Real init_step,
	func_1d & func_eval
)
{
	// Finds iterate that satisfies the strong Wolfe criteria.

	Real const param_c1( 1.0e-1 );
	Real const param_c2( 0.8 );

	Real func_value0( _func_to_beat );
	Real func_value1( 0.0 );
	Real func_value_prev( _func_to_beat );
	Real func_value_return( 0.0 );

	// Real func_from_zoom( 0.0 );

	// We already calculated the derivative to get the search direction, there's
	// no need to do it again
	Real deriv0( _deriv_sum );
	Real deriv1( 0.0 );
	Real deriv_prev( deriv0 );

	Real alpha0( 0.0 );
	Real alpha_max( 2.0 );
	Real alpha_prev( alpha0 );
	Real alpha1( init_step );
	Size iterations( 1 );

	while ( true ) {
		func_value1 = func_eval( alpha1 );

		if ( ( func_value1 > ( func_value0 + ( param_c1 * alpha1 * deriv0 ) ) ) ||
				( ( iterations > 1 ) && func_value1 >= func_value_prev ) ) {
			// Call zoom
			_last_accepted_step = zoom( alpha_prev, func_value_prev, deriv_prev, alpha1, func_value1, deriv1,
				func_value0, deriv0, func_value_return, func_eval );
			store_current_derivatives( func_eval._dE_dvars );
			return func_value_return;
		}

		deriv1 = func_eval.dfunc( alpha1 );

		if ( std::abs( deriv1 ) <= -1.0 * param_c2 * deriv0 ) {
			//  if ( deriv1 > param_c2 * deriv0 ) {
			// This corresponds to a point that satisfies both of the strong Wolfe criteria
			_last_accepted_step = alpha1;
			store_current_derivatives( func_eval._dE_dvars );
			return func_value1;
		}

		if ( deriv1 >= 0.0 ) {
			// Call zoom
			//   std::cout << "Zoom entry switch" << std::endl;
			//   std::cout << "Arguments " << alpha1 << " , " << func_value1 << " , " << deriv1 << "\n"
			//                << alpha_prev << " , " << func_value_prev << " , " << deriv_prev << "\n"
			//                << func_value_return << std::endl;
			_last_accepted_step = zoom( alpha1, func_value1, deriv1, alpha_prev, func_value_prev, deriv_prev,
				func_value0, deriv0, func_value_return, func_eval );
			store_current_derivatives( func_eval._dE_dvars );
			return func_value_return;
		}

		iterations++;
		alpha_prev = alpha1;
		func_value_prev = func_value1;
		deriv_prev = deriv1;
		alpha1 = 0.5*( alpha1 + alpha_max );
	}
}

Real
StrongWolfeLineMinimization::zoom(
	Real alpha_low,
	Real func_low,
	Real deriv_low,
	Real alpha_high,
	Real func_high,
	Real deriv_high,
	Real func_zero,
	Real deriv_zero,
	Real & func_return,
	func_1d & func_eval
)
{
	Real const param_c1( 1.0e-1 );
	Real const param_c2( 0.8 );

	Real alpha_test( 0.0 );
	Real func_test( 0.0 );
	Real deriv_test( 0.0 );

	int static step_count( 0 );
	Size iterations( 0 );
	Size const max_iterations_to_try( 8 );
	Real const min_interval_threshold( 1.0e-5 );
	Real const min_step_size( 1.0e-5 );
	Real const scale_factor( 0.66 );


	// Initial guess at test point.  Note that initially, deriv_high may be bogus, so this
	// is the only interpolation we can safely use.  Don't be tempted to try a cubic or quadratic
	// with derivative interpolations.
	alpha_test =  quadratic_interpolation( alpha_low, func_low, deriv_low, alpha_high, func_high );

	while ( true ) {

		iterations++;

		//  if( deriv_low*(alpha_high - alpha_low) > 0.0 ) {
		//   std::cout << "Bad deriv at alpha_low!" << std::endl;
		//  }

		//  std::cout << "Got point alpha_test of " << alpha_test << " from low " << alpha_low << " and high " << alpha_high << " on step " << step_count << std::endl;

		func_test = func_eval( alpha_test );
		deriv_test = func_eval.dfunc( alpha_test );

		if ( iterations > max_iterations_to_try ) {
			//   std::cout << "Warning - More'-Thuente line search exceeded max_iterations_to_try - returning current value" << std::endl;
			//   if( func_test < func_low ) {
			//    std::cout << "But at least func is less than start!" << std::endl;
			//   } else {
			//    std::cout << "And func is more than start!" << std::endl;
			//   }
			func_return = func_test;
			return alpha_test;
		}

		if ( std::abs( alpha_low - alpha_high ) <= min_interval_threshold ) {
			//   std::cout << "Warning - More'-Thuente line search interval below threshold - returning current value" << std::endl;
			//   if( func_test < func_low ) {
			//    std::cout << "But at least func is less than start!" << std::endl;
			//   } else {
			//    std::cout << "And func is more than start!" << std::endl;
			//   }
			func_return = func_test;
			return alpha_test;
		}

		if ( alpha_test <= min_step_size ) {
			//   std::cout << "Warning - More'-Thuente line search interval below threshold - returning current value" << std::endl;
			//   if( func_test < func_low ) {
			//    std::cout << "But at least func is less than start!" << std::endl;
			//   } else {
			//    std::cout << "And func is more than start!" << std::endl;
			//   }
			func_return = func_test;
			return alpha_test;
		}

		step_count++;

		if ( ( func_test > func_zero + param_c1 * alpha_test * deriv_zero ) ||
				( func_test >= func_low ) ) {

			// This is the case where the value at the test point is too high,
			// so we need to find a new test point in between the start point and
			// this point via interpolation.

			//   std::cout << "Branch/Case 1" << std::endl;

			alpha_high = alpha_test;
			func_high = func_test;
			deriv_high = deriv_test;

			Real cubic_interp = cubic_interpolation( alpha_low, func_low, deriv_low, alpha_test, func_test, deriv_test );
			Real quadratic_interp = quadratic_interpolation( alpha_low, func_low, deriv_low, alpha_test, func_test );

			if ( std::abs( cubic_interp - alpha_low ) < std::abs( quadratic_interp - alpha_low ) ) {
				alpha_test = cubic_interp;
			} else {
				alpha_test = 0.5 * ( cubic_interp + quadratic_interp );
			}

		} else {
			if ( std::abs( deriv_test ) <= -1.0 * param_c2 * deriv_zero ) {
				//   if( deriv_test > param_c2 * deriv_zero ) {
				// If we hit this we are done
				func_return = func_test;
				return alpha_test;
			}

			Real save_alpha_test = alpha_test;
			Real save_func_test = func_test;
			Real save_deriv_test = deriv_test;


			// Logic for next test point - in each case the value at the test point was lower than
			// the start point, but the Wolfe curvature criteria was not met.
			if ( deriv_test*deriv_low < 0.0 ) {

				// This case is straightforward.  The derivative at the test point is opposite to that
				// of the start point, so the bracketing is good.  We calculate the cubic and quad w/ deriv
				// interpolated guesses at the minimum, then take whichever is closest to the test point.
				// We want to be closer to the test point since it satisfies the sufficient descent
				// criterion, and we just need to get closer to the actual minimum to decrease the
				// derivative an satisfy the curvature criterion.

				//    std::cout << "Case 2" << std::endl;

				Real cubic_interp = cubic_interpolation( alpha_low, func_low, deriv_low, alpha_test, func_test, deriv_test );
				Real quadratic_interp = secant_interpolation( alpha_low, deriv_low, alpha_test, deriv_test );

				//    std::cout << "Got cubic alpha_test of " << cubic_interp << " from low " << alpha_low << " and test " << alpha_test << " on step " << step_count << std::endl;
				//    std::cout << "Got secant alpha_test of " << quadratic_interp << " from low " << alpha_low << " and test " << alpha_test << " on step " << step_count << std::endl;

				//    std::cout << "Derivs are low " << deriv_low << " and test " << deriv_test << std::endl;

				if ( std::abs( cubic_interp - alpha_test ) >= std::abs( quadratic_interp - alpha_test ) ) {
					//     std::cout << "Using cubic" << std::endl;
					alpha_test = cubic_interp;
				} else {
					//     std::cout << "Using quadratic" << std::endl;
					alpha_test = quadratic_interp;
				}

			} else if ( std::abs( deriv_test ) <= std::abs( deriv_low ) ) {
				// This is the hardest case - the slopes at low and test are in the same direction
				// and func_test is lower than func_low, so we're going in the correct direction, but
				// it's hard to tell how far to go next.

				// This is not as rigorous as in More' and Thuente 1994

				//    std::cout << "Case 3" << std::endl;

				//    std::cout << "Low : " << alpha_low << " , " << func_low << " , " << deriv_low << std::endl;
				//    std::cout << "Test: " << alpha_test << " , " << func_test << " , " << deriv_test << std::endl;
				//    std::cout << "High: " << alpha_high << " , " << func_high << " , " << deriv_high << std::endl;

				Real try_alpha_test = secant_interpolation( alpha_low, deriv_low, alpha_test, deriv_test );
				//    Real try_alpha_test = quadratic_interpolation( alpha_test, func_test, deriv_test, alpha_high, func_high );
				// If this is outside of the bounds of (alpha_test, alpha_high), scale back
				Real scaled_alpha_test = alpha_test + scale_factor*( alpha_high - alpha_test );

				//    std::cout << "Got interpolated alpha_test of " << try_alpha_test << std::endl;
				if ( alpha_test > alpha_low ) {
					//     std::cout << "Rescaling case 1!" << std::endl;
					alpha_test = std::min( scaled_alpha_test, try_alpha_test );
				} else {
					//     std::cout << "Rescaling case 2!" << std::endl;
					alpha_test = std::max( scaled_alpha_test, try_alpha_test );
				}
				//    std::cout << "Post out-of-bounds check alpha_test is " << alpha_test << std::endl;

			} else {
				//    std::cout << "Case 4" << std::endl;

				// deriv_high may be totally bogus?
				if ( deriv_high == 0.0 ) {
					alpha_test = quadratic_interpolation( alpha_test, func_test, deriv_test, alpha_high, func_high );
				} else {
					alpha_test = cubic_interpolation( alpha_test, func_test, deriv_test, alpha_high, func_high, deriv_high );
				}
			}


			// Finally, adjust boundaries of the interval

			if ( ( deriv_test * ( alpha_test - alpha_low ) ) >= 0.0 ) {
				//    std::cout << "Branch 2 if check" << std::endl;
				alpha_high = alpha_low;
				func_high = func_low;
				deriv_high = deriv_low;
			}
			//   std::cout << "Branch 2" << std::endl;
			alpha_low = save_alpha_test;
			func_low = save_func_test;
			deriv_low = save_deriv_test;
		}

		// Safeguard the new step
		if ( alpha_high > alpha_low ) {
			alpha_test = std::min( alpha_low + scale_factor*( alpha_high - alpha_low ), alpha_test );
		} else {
			alpha_test = std::max( alpha_low + scale_factor*( alpha_high - alpha_low ), alpha_test );
		}

		//  std::cout << "For next round" << std::endl;
		//  std::cout << "Low : " << alpha_low << " , " << func_low << " , " << deriv_low << std::endl;
		//  std::cout << "Test: " << alpha_test << " , " << func_test << " , " << deriv_test << std::endl;
		//  std::cout << "High: " << alpha_high << " , " << func_high << " , " << deriv_high << std::endl;

	}

}

void
LineMinimizationAlgorithm::store_current_derivatives(
	Multivec & curr_derivs
)
{
	// std::cout << "Storing derivatives" << " size " << _stored_derivatives.size() << std::endl;
	for ( uint i =  1 ; i <= _stored_derivatives.size() ; ++i ) {
		_stored_derivatives[i] = curr_derivs[i];
		//  if( i <= 10 ) std::cout << "Storing deriv " << i << " as " << _stored_derivatives[i] << std::endl;
	}
}

void
LineMinimizationAlgorithm::fetch_stored_derivatives(
	Multivec & set_derivs
)
{
	// std::cout << "Fetching derivatives" << " size " << _stored_derivatives.size() << std::endl;
	for ( uint i =  1 ; i <= _stored_derivatives.size() ; ++i ) {
		set_derivs[i] = _stored_derivatives[i];
		//  if( i <= 10 ) std::cout << "Fetching deriv " << i << " as " << set_derivs[i] << std::endl;
	}
}


Real
LineMinimizationAlgorithm::quadratic_interpolation(
	Real alpha_low,
	Real func_low,
	Real deriv_low,
	Real alpha_high,
	Real func_high
)
{
	return ( 2.0*alpha_low*( func_high - func_low ) -
		deriv_low * ( alpha_high*alpha_high - alpha_low*alpha_low ) ) /
		( 2.0 * ( (func_high - func_low ) - ( deriv_low * ( alpha_high - alpha_low ) ) ) );
}

Real
LineMinimizationAlgorithm::quadratic_deriv_interpolation(
	Real alpha_low,
	Real func_low,
	Real deriv_low,
	Real alpha_high,
	Real func_high,
	Real // deriv_high
)
{
	// Real const alpha_diff( alpha_high - alpha_low );
	// Real const deriv_diff( deriv_high - deriv_low );
	// return ( alpha_low  - ( deriv_low * alpha_diff / deriv_diff ) );
	return ( alpha_low + ( (deriv_low/((func_low - func_high)/(alpha_high - alpha_low) + deriv_low))*0.5) *(alpha_high -alpha_low) );
}

Real
LineMinimizationAlgorithm::secant_interpolation(
	Real alpha_low,
	Real deriv_low,
	Real alpha_high,
	Real deriv_high
)
{
	// Real const alpha_diff( alpha_high - alpha_low );
	// Real const deriv_diff( deriv_high - deriv_low );
	return ( alpha_low + (deriv_low/(deriv_low - deriv_high))*(alpha_high -alpha_low) );
}

Real
LineMinimizationAlgorithm::cubic_interpolation(
	Real alpha_low,
	Real func_low,
	Real deriv_low,
	Real alpha_high,
	Real func_high,
	Real deriv_high
)
{
	Real cubic_beta1 = deriv_low + deriv_high - 3.0*( func_low - func_high )/(alpha_low - alpha_high);
	Real max_beta_derivs = std::max( std::max( std::abs( cubic_beta1), std::abs( deriv_low ) ), std::abs( deriv_high ) );
	Real gamma = max_beta_derivs * std::sqrt( std::pow( cubic_beta1 / max_beta_derivs, 2.0 ) - (deriv_low/max_beta_derivs)*(deriv_high/max_beta_derivs));
	if ( alpha_high < alpha_low ) gamma = -gamma;
	Real numer = (gamma - deriv_low) + cubic_beta1;
	Real denom = ((gamma - deriv_low) + gamma) + deriv_high;
	Real fraction = numer/denom;
	// Real cubic_beta2 = std::sqrt( cubic_beta1*cubic_beta1 - (deriv_low * deriv_high) );
	// return  ( alpha_high - (alpha_high - alpha_low)*( deriv_high + cubic_beta2 - cubic_beta1 ) /
	//              (deriv_high - deriv_low + 2.0*cubic_beta2 ) );
	return  ( alpha_low + fraction * ( alpha_high - alpha_low ) );
}


} // namespace optimization
} // namespace core
