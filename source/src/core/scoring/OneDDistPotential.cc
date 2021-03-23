// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/OneDDistPotential.cc
/// @brief One-dimensional distance potential
/// @detailed
/// @author Andy Watkins (amw579stanford.edu)


#include <core/scoring/OneDDistPotential.hh>
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
#include <utility/numbers.hh> // for isfinite

static basic::Tracer TR( "core.scoring.OneDDistPotential" );

using namespace utility;
using namespace utility::tools;
using namespace utility::json_spirit;
using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace core {
namespace scoring {

//Constructor
OneDDistPotential::OneDDistPotential( std::string const & filename ):
	boundary_( numeric::interpolation::FLAT )
{
	mObject json;
	read_tensor_from_file( filename, tensor_, json );
	initialize_from_json( json );
}

OneDDistPotential::OneDDistPotential( numeric::MathNTensor< core::Real, 1 > const & tensor,
	utility::json_spirit::mObject const & json ):
	boundary_( numeric::interpolation::FLAT )
{
	tensor_ = tensor;
	initialize_from_json( json );
}

//Destructor
OneDDistPotential::~OneDDistPotential() = default;

void
inline
get_minval_binwidth( numeric::MathNTensor< core::Real, 1 > const & T,
	utility::json_spirit::mObject const & json,
	utility::fixedsizearray1< numeric::Real, 1 > & minval,
	utility::fixedsizearray1< numeric::Real, 1 > & binwidth)
{
	using namespace utility::json_spirit;
	mArray json_binwidth( get_mArray( json, "binwidth" ) );
	mArray json_minval( get_mArray( json, "minval" ) );
	mArray json_maxval( get_mArray( json, "maxval" ) );

	// TR << "For this json..." << std::endl;
	// TR << "Max val was " << json_maxval[0].get_real() << std::endl;
	// TR << "Min val was " << json_minval[0].get_real() << std::endl;
	// TR << "Bin wid was " << json_binwidth[0].get_real() << std::endl;
	// TR << "Ratio is " << numeric::Size( ( json_maxval[0].get_real() - json_minval[0].get_real() )/json_binwidth[0].get_real() + 1 ) << std::endl;
	// TR << "nbins is " << T.n_bins( 1 ) << std::endl;

	runtime_assert( numeric::Size( ( json_maxval[0].get_real() - json_minval[0].get_real() )/json_binwidth[0].get_real() + 1 ) == T.n_bins( 1 ) );
	minval[ 1 ]   = json_minval[0].get_real();
	binwidth[ 1 ] = json_binwidth[0].get_real();

}

///////////////////////////////////////////////////////////////////////////////////
void
OneDDistPotential::initialize_from_json( utility::json_spirit::mObject const & json ) {
	get_minval_binwidth( tensor_, json, minval_, binwidth_ );
	// Rmax_ = get_real_or_zero( json, "Rmax" );
	// Emax_ = get_real_or_zero( json, "Emax" );
}

///////////////////////////////////////////////////////////////////////////////////
core::Real
OneDDistPotential::evaluate( core::Real const dist ) const
{
	core::Real dummy_deriv;
	return evaluate( dist, false, dummy_deriv );
}

///////////////////////////////////////////////////////////////////////////////////
core::Real
OneDDistPotential::get_derivative( core::Real const dist ) const
{
	core::Real deriv;
	evaluate( dist, true, deriv );
	return deriv;
}

core::Real
OneDDistPotential::evaluate( core::Real const dist,
	bool const compute_deriv, core::Real & deriv ) const
{
	// Real const R( t.length() );
	// if ( Rmax_ > 0.0 && R > Rmax_ ) {
	//  if ( compute_deriv ) deriv = evaluate_constraining_potential_derivative( t );
	//  return evaluate_constraining_potential( R );
	// }

	utility::fixedsizearray1< core::Real, 1 > outvals( 0.0 ), tensor_deriv( 0.0 );
	outvals[ 1 ] = dist;

	Real value( 0.0 );
	if ( dist > minval_[ 1 ] && dist < minval_[ 1 ] + binwidth_[ 1 ] * ( tensor_.n_bins()[ 1 ] - 3 )  ) {
		value = numeric::interpolation::polycubic_interpolate_catmull_rom( tensor_, minval_, binwidth_, outvals, boundary_, tensor_deriv, compute_deriv );
		if ( value < -6 ) {
			TR << "dist " << dist << " " << " minval " << minval_[ 1 ] << " " <<  binwidth_[ 1 ] << " nbins " << tensor_.n_bins()[ 1 ] << std::endl;
		}
	}
	runtime_assert( utility::isfinite( value ) );

	if ( compute_deriv ) {
		deriv = tensor_deriv[1];
	}
	return value;
}

////////////////////////////////////////////////////////////////////////////////////
/// @details super-steep potential beyond Rmax to really penalize movement outside Rmax.
///          Arbitrarily using R^6 potential.
// core::Real
// OneDDistPotential::evaluate_constraining_potential( Distance const & R ) const
// {
//  return ( Emax_ * std::pow( ( R / Rmax_ ), 6.0 ) );
// }

// ////////////////////////////////////////////////////////////////////////////////////
// std::pair< Vector, Vector >
// OneDDistPotential::evaluate_constraining_potential_derivative( Vector const & t ) const
// {
//  Distance const R( t.length() );
//  Vector const r( t / R );
//  Vector const trans_deriv ( r * Emax_ * ( 6.0 / Rmax_) * std::pow( ( R / Rmax_ ), 5.0 ) );
//  Vector const rot_deriv( 0.0 );
//  return std::make_pair( trans_deriv, rot_deriv );
// }

} //scoring
} //core
