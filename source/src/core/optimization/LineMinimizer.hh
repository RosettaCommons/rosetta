// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/LineMinimizer.hh
/// @brief  line minimizer classes
/// @author Phil Bradley
/// @author Jim Havranek


#ifndef INCLUDED_core_optimization_LineMinimizer_hh
#define INCLUDED_core_optimization_LineMinimizer_hh


// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>

#include <core/optimization/LineMinimizer.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

namespace core {
namespace optimization {

// Functor that stores current position and
// search direction to provide a Real->Real mapping
// for univariate minimization

class func_1d {
private:
	Multivec const _starting_point;
	Multivec const _search_direction;
	Multivec _eval_point;
public:
	Multivec _dE_dvars;
private:
	Multifunc const & _func;
	int _eval_count;
	int _deriv_count;
	// Real _search_direction_magnitude;
public:
	func_1d( Multivec & start, Multivec & dir, Multifunc const & score_fxn ) :
		_starting_point( start ),
		_search_direction( dir ),
		_eval_point( dir.size(), 0.0 ),
		_dE_dvars( dir.size(), 0.0 ),
		_func( score_fxn ),
		_eval_count( 0 ),
		_deriv_count( 0 )
		// _search_direction_magnitude( 0.0 )
	{
		debug_assert( _starting_point.size() == _search_direction.size() );
		//for( uint i =  1 ; i <= _starting_point.size() ; ++i ) {
		// _search_direction_magnitude += _search_direction[i] * _search_direction[i];
		//}
		//_search_direction_magnitude = std::sqrt( _search_direction_magnitude );
	};

	Real operator() ( Real displacement ) {
		_eval_count++;
		for ( uint i =  1 ; i <= _starting_point.size() ; ++i ) {
			_eval_point[i] = _starting_point[i] +
				( displacement * _search_direction[i] );
		}
		return _func( _eval_point );
	};

	Real dfunc( Real displacement ) {
		_deriv_count++;
		Real dot_product( 0.0 );
		for ( uint i =  1 ; i <= _starting_point.size() ; ++i ) {
			_eval_point[i] = _starting_point[i] +
				( displacement * _search_direction[i] );
		}
		_func.dfunc( _eval_point, _dE_dvars );
		for ( uint i =  1 ; i <= _starting_point.size() ; ++i ) {
			dot_product += ( _dE_dvars[i] * _search_direction[i] );
		}
		return dot_product;
	};

	void reset_eval_count() { _eval_count = 0; };
	int get_eval_count() { return _eval_count; };
	int get_deriv_count() { return _deriv_count; };
	/// @brief Error condition wherein the computed gradient does not match the actual gradient;
	/// invokes the Multifunc::dump( vars, vars2 ) method.
	void dump( Real displacement );
	//Real search_direction_magnitude() { return _search_direction_magnitude; };
};

/////////////////////////////////////////////////////////
// base class / interface for line minimizers
/////////////////////////////////////////////////////////
class LineMinimizationAlgorithm : public utility::pointer::ReferenceCount
{
public:
	~LineMinimizationAlgorithm() override;
	LineMinimizationAlgorithm( Multifunc const & score_fxn, Size dimension ) :
		_stored_derivatives( dimension, 0.0 ),
		_last_accepted_step( 1.0 ), _func_to_beat( 0.0 ),
		_deriv_sum( 0.0 ),_num_linemin_calls(0),_tolerance( 0.1 ),
		_func( score_fxn ), _nonmonotone( false ), _silent( false ) {};
	virtual Real operator()( Multivec & , Multivec & ){ return 0.0; };
	virtual bool provide_stored_derivatives(){ return false; };
	bool nonmonotone() { return _nonmonotone; };
	void store_current_derivatives( Multivec & curr_derivs );
	void fetch_stored_derivatives( Multivec & get_derivs );
	Real quadratic_interpolation( Real point1, Real func1, Real deriv1, Real point2, Real func2 );
	Real quadratic_deriv_interpolation( Real point1, Real func1, Real deriv1, Real point2, Real func2, Real deriv2 );
	Real secant_interpolation( Real point1, Real deriv1, Real point2, Real deriv2 );
	Real cubic_interpolation( Real point1, Real func1, Real deriv1, Real point2, Real func2, Real deriv2 );

	bool silent() { return _silent; };
	void silent(bool s_in) { _silent=s_in; };


	Multivec _stored_derivatives;
	Real _last_accepted_step;
	Real _func_to_beat;
	Real _deriv_sum;
	int _num_linemin_calls;
protected:
	Real const _tolerance;
	Multifunc const & _func;
	bool _nonmonotone;
	bool _silent;
};

/////////////////////////////////////////////////////////
// concrete line minimizer - Brent's method
/////////////////////////////////////////////////////////

class BrentLineMinimization : public LineMinimizationAlgorithm
{
public:
	BrentLineMinimization( Multifunc const & score_fxn, Size dim ) :
		LineMinimizationAlgorithm( score_fxn, dim ),
		_ax( 0.0 ), _bx( 0.2 ), _xx( 0.1 ), _abs_tolerance( 0.01 ),
		deriv_cutoff_(0.0001)
	{};
	Real operator()( Multivec & curr_pos, Multivec & curr_dir ) override;
	void set_deriv_cutoff( core::Real const &val ) { deriv_cutoff_=val; }
	void MNBRAK( Real & AX, Real & BX, Real & CX, Real & FA, Real & FB,
		Real & FC, func_1d & func_eval) const;
	Real BRENT( Real const AX, Real const BX, Real const CX, Real & FA,
		Real & FB, Real const FC, Real const TOL, func_1d & func_eval);
	Real _ax;
	Real _bx;
	Real _xx;
	Real _abs_tolerance;
	Real deriv_cutoff_;
};

/////////////////////////////////////////////////////////
// concrete line minimizer - Armijo's method
/////////////////////////////////////////////////////////

class ArmijoLineMinimization : public LineMinimizationAlgorithm
{
public:
	ArmijoLineMinimization( Multifunc const & score_fxn, bool nonmonotone, Size dim ) :
		LineMinimizationAlgorithm( score_fxn, dim ),
		_num_calls( 0 ) { _nonmonotone = nonmonotone; };
	bool provide_stored_derivatives() override{ return false; };
	Real operator()( Multivec & curr_pos, Multivec & curr_dir ) override;
	Real Armijo( Real init_step, func_1d & func_eval );

	int _num_calls;
};

/////////////////////////////////////////////////////////
// concrete line minimizer - Satisfies strong Wolfe conditions
// Roughly following More' amd Thuente, 1994
/////////////////////////////////////////////////////////

class StrongWolfeLineMinimization : public LineMinimizationAlgorithm
{
public:
	StrongWolfeLineMinimization( Multifunc const & score_fxn, bool nonmonotone, Size dim ) :
		LineMinimizationAlgorithm( score_fxn, dim ),
		_nonmonotone( nonmonotone ),
		_num_calls( 0 ) {};
	bool provide_stored_derivatives() override{ return true; };
	Real operator()( Multivec & curr_pos, Multivec & curr_dir ) override;
	Real StrongWolfe( Real init_step, func_1d & func_eval );
	Real zoom( Real alpha_low, Real func_low, Real deriv_low, Real alpha_high, Real func_high, Real deriv_high,
		Real func_zero, Real deriv_zero, Real & func_return, func_1d & func_eval );
	bool _nonmonotone;
	int _num_calls;
};

} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_LineMinimizer_HH
