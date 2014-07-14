// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/Minimizer.hh
/// @brief  Simple low-level minimizer class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_Minimizer_hh
#define INCLUDED_core_optimization_Minimizer_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/Multifunc.hh>
// AUTO-REMOVED #include <core/optimization/LineMinimizer.hh>

#include <core/optimization/LineMinimizer.fwd.hh>
#include <core/optimization/Multifunc.fwd.hh>
#include <utility/vector1.hh>



namespace core {
namespace optimization {


/**************************************************
 *
 * Rough outline of how to structure this:
 *
 * Make a base 'minimizer' class
 *
 * Sub-class into univariate and multivariate minimizers
 * -> actually, could we treat linmin as an instance of multivariate
 *  	minimization, with a single pass of steepest descent
 *
 * The trick is how to mix and match convergence criteria, descent
 * direction generation, and line minimization schemes
 *
 * convergence criteria could be a function or a functor.  Descent
 * direction algorithms probably need to be functors, since they
 * have different storage needs.
 *
 *
 *
 **********************************************/

//**************************************
//*** Begin convergence test section ***
//**************************************

// base convergence test class
class ConvergenceTest {
public:
	virtual bool operator()( Real Fnew, Real Fold ) = 0;
	virtual ~ConvergenceTest() {}
};

// concrete convergence test class - classic dfpmin
class DFPMinConvergedFractional : public ConvergenceTest {
public:
	DFPMinConvergedFractional( Real _tol, Real _eps = 1.0e-10 ) : tolerance( _tol ), eps( _eps ) {};
	virtual ~DFPMinConvergedFractional() {}
	virtual bool operator()( Real Fnew, Real Fold );
private:
	Real tolerance;
	Real eps;
};

// concrete convergence test class - "atol" dfpmin
class DFPMinConvergedAbsolute : public ConvergenceTest {
public:
	DFPMinConvergedAbsolute( Real _tol ) : tolerance( _tol ) {}
	virtual ~DFPMinConvergedAbsolute() {}
	virtual bool operator()( Real Fnew, Real Fold );
private:
	Real tolerance;
};

//**************************************
//***  End  convergence test section ***
//**************************************

//***************************************
//*** Begin descent direction section ***
//***************************************

// base convergence test class
class DescentDirectionAlgorithm {
public:
	DescentDirectionAlgorithm( Multifunc const & ) {} // func_ is unused? what's this class for? : func_( in_func_ ) {}
//	Multivec operator()(){};
	void initialize(){};
private:
	// Multifunc const & func_;
};

//********************************************
//*** Begin multivariate minimizer section ***
//********************************************

class JJH_Minimizer {
public:
	JJH_Minimizer(
		Multifunc const & score_fxn,
		LineMinimizationAlgorithm & line_min_alg,
		ConvergenceTest & converge_test,
		DescentDirectionAlgorithm & desc_dir ):
		_func( score_fxn ), _line_min( line_min_alg ),
		_converged( converge_test ), _get_direction( desc_dir ){}
	Real run( Multivec & angles );
private:
	Multifunc const & _func;
	LineMinimizationAlgorithm & _line_min;
	ConvergenceTest & _converged;
	DescentDirectionAlgorithm &	_get_direction;
};

//********************************************
//***  End  multivariate minimizer section ***
//********************************************

//
// Data stored at each LBFGS iteration
class lbfgs_iteration_data {
public:
	core::Real alpha;
	Multivec s;
	Multivec y;
	core::Real ys; // = dot(y,s)
};


/// @brief Simple low-level minimizer class
class Minimizer {

public:
	Minimizer( Multifunc & func_in, MinimizerOptions const & options_in );

	Real
	run( Multivec & phipsi_inout );

private:
	void
	linmin(
		Multivec & P,
		Multivec & XI,
		Real & FRET,
		int const ITMAX
	) const;

	void
	dfpmin(
		Multivec & P,
		Real & FRET,
		ConvergenceTest & converge_test,
		int const ITMAX
	) const;

	void
	dfpmin_armijo(
		Multivec & P,
		Real & FRET,
		ConvergenceTest & converge_test,
		LineMinimizationAlgorithmOP line_min,
		int const ITMAX
	) const;

	void
	lbfgs(
		Multivec & P,
		Real & FRET,
		ConvergenceTest & converge_test,
		LineMinimizationAlgorithmOP line_min,
		int const ITMAX,
		bool w_rescore = false
	) const;

	Multifunc & func_;
	MinimizerOptions options_;
}; // Minimizer

} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_Minimizer_HH
