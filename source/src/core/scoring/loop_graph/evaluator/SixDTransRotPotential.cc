// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/SixDTransRotPotential.cc
/// @brief Six-dimensional potential (translations, rotations) that can be evaluated on a jump
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.hh>
#include <core/scoring/loop_graph/evaluator/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/MathNTensor_io.hh>
#include <numeric/interpolation/interpolation.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

static basic::Tracer TR( "core.scoring.loop_graph.evaluator.SixDTransRotPotential" );

using namespace utility;
using namespace utility::tools;
using namespace utility::json_spirit;
using namespace basic::options;
using namespace basic::options::OptionKeys;

///////////////////////////////////////////////////////////////////////////////////
// @detailed
//
// Six dimensions are:
//
//   x, y, z    --> translation vector, in Angstroms
//   vx, vy, vz --> rotation vector, in degrees
//
// vx, vy, and vz define axis and angle [e.g., they are the 'exponential map' of SO(3)].
// Defined for the sphere |v| < 180.0 degrees  (i.e., boundary of pi radians).
//
// The potentials can be defined through, e.g., Monte Carlo calculations. Note that,
//   when using rotation vector representations, the 'reference' distribution is not
//   uniform in vx, vy, vx but instead is weighted like:
//
//     [ sin( |v|/2 ) / (|v|/2) ]^2
//
//  (this equation has v in radians)
//
// This integrates to 8 pi^2, where |v| (the magnitude of the rotation) fills out the sphere
//   from 0 out to pi radians (180 degrees). There are derivations of this
//   volume element throughout the literature, e.g., R.E. Miles, Biometrika
//  Vol. 52, No. 3/4 (Dec., 1965), pp. 636-639. The simplest derivation I've come up
//  with is written out notes posted in github below.
//
// Added option to enforce continuity at |v| = 180.0, with v wrapping around to
//  v - 2 pi v.normalized().
//
// Originally used polylinear interpolation because it takes a huge amount of memory
//  to precalculate 6D polycubic spline; now using 'on-the-fly' Catmull-Rom cubic interpolation by default.
//
// More information at:
//
//  https://github.com/rhiju/loop_close/blob/master/notes/Das_6DPotential_Derivs.pdf
//  https://github.com/rhiju/loop_close/blob/master/notes/Das_VolumeElementSO3_AxisAngle.pdf
//
// and command-lines are in that repo too. See also,
//
//  "Quaternions in molecular modeling", Karney, J Mol. Graph. Mod. 25 (2007) 595-604
//
///////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

//Constructor
SixDTransRotPotential::SixDTransRotPotential( std::string const & filename ):
	enforce_continuity_at_pi_( true ),
	turn_off_rotation_dependence_( false ),
	use_cubic_interp_( option[ score::loop_close::use_cubic_interp ]() ),
	boundary_( numeric::interpolation::FLAT )
{
	mObject json;
	read_tensor_from_file( filename, tensor_, json );
	initialize_from_json( json );
}

SixDTransRotPotential::SixDTransRotPotential( numeric::MathNTensor< core::Real, 6 > const & tensor,
	utility::json_spirit::mObject const & json ):
	enforce_continuity_at_pi_( true ),
	turn_off_rotation_dependence_( false ),
	use_cubic_interp_( option[ score::loop_close::use_cubic_interp ]() ),
	boundary_( numeric::interpolation::FLAT )
{
	tensor_ = tensor;
	initialize_from_json( json );
}

//Destructor
SixDTransRotPotential::~SixDTransRotPotential()
{}

///////////////////////////////////////////////////////////////////////////////////
void
SixDTransRotPotential::initialize_from_json( utility::json_spirit::mObject const & json ) {
	get_minval_binwidth( tensor_, json, minval_, binwidth_ );
	Rmax_ = get_real_or_zero( json, "Rmax" );
	Emax_ = get_real_or_zero( json, "Emax" );
}

///////////////////////////////////////////////////////////////////////////////////
core::Real
SixDTransRotPotential::evaluate( core::kinematics::Jump const & j ) const
{
	Vector const & t( j.get_translation() );
	Vector const rotation_vector( numeric::rotation_axis_angle( j.get_rotation() ) * (180.0 / numeric::constants::r::pi) );
	std::pair< Vector, Vector > dummy_deriv;
	return evaluate( t, rotation_vector, false, dummy_deriv );
}

///////////////////////////////////////////////////////////////////////////////////
std::pair< Vector, Vector >
SixDTransRotPotential::get_derivative( core::kinematics::Jump const & j ) const
{
	Vector const & t( j.get_translation() );
	Vector const rotation_vector( numeric::rotation_axis_angle( j.get_rotation() ) * (180.0 / numeric::constants::r::pi) );
	std::pair< Vector, Vector > deriv;
	evaluate( t, rotation_vector, true, deriv );
	return deriv;
}

void
zero_out_rotation_components( utility::fixedsizearray1< core::Real, 6 > & x ) {
	x[ 4 ] = 0.0;
	x[ 5 ] = 0.0;
	x[ 6 ] = 0.0;
}

///////////////////////////////////////////////////////////////////////////////////
/// @details magnitude of rotation_vector should be rotation angle in degrees, maximum of 180.
///          derivative is pair corresponding to translation derivative and rotvector derivative
core::Real
SixDTransRotPotential::evaluate( Vector const & t, Vector const & rotation_vector,
	bool const compute_deriv, std::pair< Vector, Vector > & deriv ) const
{
	Real const R( t.length() );
	if ( Rmax_ > 0.0 && R > Rmax_ ) {
		if ( compute_deriv ) deriv = evaluate_constraining_potential_derivative( t );
		return evaluate_constraining_potential( R );
	}

	utility::fixedsizearray1< core::Real, 6 > outvals( 0.0 ), tensor_deriv( 0.0 );
	outvals[ 1 ] = t.x();
	outvals[ 2 ] = t.y();
	outvals[ 3 ] = t.z();
	outvals[ 4 ] = rotation_vector.x();
	outvals[ 5 ] = rotation_vector.y();
	outvals[ 6 ] = rotation_vector.z();

	if ( turn_off_rotation_dependence_ ) zero_out_rotation_components( outvals ); // for checking derivatives.

	Real value( 0.0 );
	if ( use_cubic_interp_ ) {
		value = numeric::interpolation::polycubic_interpolate_catmull_rom( tensor_, minval_, binwidth_, outvals, boundary_, tensor_deriv, compute_deriv );
	} else {
		value = numeric::interpolation::multilinear_interpolation( tensor_, minval_, binwidth_, outvals, tensor_deriv, compute_deriv );
	}
	if ( enforce_continuity_at_pi_ ) adjust_near_pi( rotation_vector, value, outvals, compute_deriv, tensor_deriv );

	if ( turn_off_rotation_dependence_ ) zero_out_rotation_components( tensor_deriv ); // for checking derivatives.

	if ( compute_deriv ) {
		deriv = std::make_pair( Vector( tensor_deriv[1], tensor_deriv[2], tensor_deriv[3] ),
			Vector( tensor_deriv[4], tensor_deriv[5], tensor_deriv[6] ) );
	}
	return value;
}

///////////////////////////////////
void
SixDTransRotPotential::adjust_near_pi( Vector const & rotation_vector,
	core::Real & value,
	utility::fixedsizearray1< core::Real, 6 > & outvals,
	bool const & compute_deriv,
	utility::fixedsizearray1< core::Real, 6 > & tensor_deriv ) const
{
	if ( std::abs( rotation_vector.length() - 180.0 ) < binwidth_[4] ) { // in the danger zone.
		Vector const rotation_vector_flipped = rotation_vector - 360.0 * rotation_vector.normalized();
		outvals[ 4 ] = rotation_vector_flipped.x();
		outvals[ 5 ] = rotation_vector_flipped.y();
		outvals[ 6 ] = rotation_vector_flipped.z();
		utility::fixedsizearray1< core::Real, 6 > tensor_deriv2;
		Real value2 = numeric::interpolation::multilinear_interpolation( tensor_, minval_, binwidth_, outvals, tensor_deriv2, compute_deriv );
		Real a = 1.0 - ( 180.0 - rotation_vector.length() ) / binwidth_[ 4 ]; // how close are we to pi? Goes from 0.0 to 1.0.
		runtime_assert( a >= 0.0 );
		if ( a - 1.0 > 1.0e-3 ) {
			std::cerr << "rotation_vector is " << rotation_vector << "  and has length " << rotation_vector.length() << "; a is: " << a << std::endl;
		}
		runtime_assert( (a - 1.0) <= 1.0e-3 );
		value = (1.0 - 0.5 * a ) * value + ( 0.5 * a ) * value2; // at pi, use half-half
		if ( compute_deriv ) {
			for ( Size n = 1; n <= 6; n++ ) tensor_deriv[n] = (1.0 - 0.5 * a ) * tensor_deriv[n] + ( 0.5 * a ) * tensor_deriv2[n];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
/// @details super-steep potential beyond Rmax to really penalize movement outside Rmax.
///          Arbitrarily using R^6 potential.
core::Real
SixDTransRotPotential::evaluate_constraining_potential( Distance const & R ) const
{
	return ( Emax_ * std::pow( ( R / Rmax_ ), 6.0 ) );
}

////////////////////////////////////////////////////////////////////////////////////
std::pair< Vector, Vector >
SixDTransRotPotential::evaluate_constraining_potential_derivative( Vector const & t ) const
{
	Distance const R( t.length() );
	Vector const r( t / R );
	Vector const trans_deriv ( r * Emax_ * ( 6.0 / Rmax_) * std::pow( ( R / Rmax_ ), 5.0 ) );
	Vector const rot_deriv( 0.0 );
	return std::make_pair( trans_deriv, rot_deriv );
}

} //evaluator
} //loop_graph
} //scoring
} //core
