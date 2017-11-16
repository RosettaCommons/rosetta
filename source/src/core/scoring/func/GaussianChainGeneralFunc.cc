// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainGeneralFunc.hh
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainGeneralFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <numeric/constants.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <iostream>
#include <cmath>
#ifdef WIN32
#include <boost/math/special_functions/erf.hpp>
#endif

// See GaussianChainFunc.cc for more information, including link to mathematical derivation.

using numeric::constants::d::pi;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.scoring.func.GaussianChainGeneralFunc" );

namespace core {
namespace scoring {
namespace func {

GaussianChainGeneralFunc::GaussianChainGeneralFunc(
	Real const gaussian_variance,
	Real const loop_fixed_cost,
	utility::vector1< Real > const & other_distances ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost ),
	other_distances_( other_distances )
{
	initialize_parameters();
}

FuncOP
GaussianChainGeneralFunc::clone() const
{
	return FuncOP( new GaussianChainGeneralFunc( gaussian_variance_, loop_fixed_cost_, other_distances_ ) );
}

bool GaussianChainGeneralFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	GaussianChainGeneralFunc const & other_downcast( static_cast< GaussianChainGeneralFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( other_distances_       != other_downcast.other_distances_       ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( loop_fixed_cost_total_ != other_downcast.loop_fixed_cost_total_ ) return false;
	if ( N_ != other_downcast.N_ ) return false;

	return true;
}

bool GaussianChainGeneralFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainGeneralFunc const * > ( &other );
}

void
GaussianChainGeneralFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	recompute_parameters();
}

void
GaussianChainGeneralFunc::recompute_parameters(){
	N_ = 1 + other_distances_.size(); // number of fixed-distance segments.

	// this is a 'prefactor' in the probability.
	loop_fixed_cost_total_ = loop_fixed_cost_;
	loop_fixed_cost_total_ += ( N_ + 2 ) * kB_T_ * log( 2.0 );
	loop_fixed_cost_total_ += kB_T_ * log( pi );
	loop_fixed_cost_total_ += ( int( N_ ) - 3 ) * -0.5 * kB_T_ * log( 2 * gaussian_variance_ );
	// further normalization factor is put into func itself.
}

//////////////////////////////////////////////
// N-th integral from x to infinity of erfc (complementary error function)
Real
GaussianChainGeneralFunc::erfc_integral( Real const & x, int const N ) const
{
	runtime_assert( N >= -3 );
	if ( N == -3 ) {
		// -3th integral (i.e., third derivative ) of erfc
		return ( 2.0 / sqrt( pi ) ) * 2 * ( -1 + 2 * x * x ) * exp( -1.0 * x * x );
	} else if ( N == -2 ) {
		// -2th integral (i.e., second derivative ) of erfc
		return ( 2.0 / sqrt( pi ) ) * 2 * x * exp( -1.0 * x * x );
	} else if ( N == -1 ) {
		// -1th integral (i.e., first derivative ) of erfc
		return ( 2.0 / sqrt( pi ) ) * exp( -1.0 * x * x );
	} else if ( N == 0 ) {
		// zeroth integral of erfc, i.e. erfc
#ifdef WIN32
		return boost::math::erfc( x );
#else
		return erfc( x );
#endif
	} else {
		// N-th integral can be derived from a recurrence if we know N-1-th and N-2-th integrals
		// (Note: easy to derive via integration by parts.)
		// Abramowitz and Stegun, p. 299
		return ( ( -x / Real( N ) ) * erfc_integral( x, N - 1 ) + ( 1 / ( 2.0 * N ) ) * erfc_integral( x, N - 2 ) );
	}
}

Real
GaussianChainGeneralFunc::func( Real const z,
	bool const calc_deriv,
	Real & deriv ) const
{
	Real const s = sqrt( 2 * gaussian_variance_ );
	Real L_sum( 0.0 ), deriv_sum( 0.0 );

	// put all distances in a vector
	utility::vector1< Distance > distances( 1, z );
	distances.append( other_distances_ );

	// Need to get all combinations of signs, e.g., if N = 4, we need the 16 multiplets, (+1,+1,+1,+1), (-1,+1,+1,+1) ...
	Size num_terms( 1 );
	for ( Size i = 1; i <= N_; i++ ) num_terms *= 2; // 2^N, e.g., 16.
	for ( Size k = 1; k <= num_terms; k++ ) {
		utility::vector1< int > sgn( N_, 0 );
		int sgn_prod( 1 );
		Real D_sum( 0.0 );
		Size count( k );
		for ( Size i = 1; i <= N_; i++ ) {
			sgn[ i ] = 2 * ( count % 2 ) - 1; // converts 0,1 to -1,1
			count /= 2; // next digit, in binary
			D_sum += sgn[ i ] * distances[ i ];
			sgn_prod *= -sgn[ i ];
		}
		D_sum /= s;
		L_sum += -1.0 * sgn_prod * erfc_integral( D_sum, int( N_ ) - 3 );
		if ( calc_deriv ) {
			// note that sign flips from L_sum because erfc_integral integrates from x to infinity, not 0 to x.
			deriv_sum += sgn[ 1 ] * sgn_prod * erfc_integral( D_sum, int( N_ ) - 4 );
		}
	}

	Real loop_energy = -kB_T_ * log( L_sum );
	for ( auto distance : distances ) loop_energy += kB_T_ * log( distance );

	if ( calc_deriv ) {
		Real const first_logderiv_term = kB_T_ / z;
		Real second_logderiv_term = deriv_sum / L_sum;
		second_logderiv_term *= -1.0 * kB_T_ / s;
		deriv = first_logderiv_term + second_logderiv_term;
	}
	return ( loop_fixed_cost_total_ + loop_energy );
}

Real
GaussianChainGeneralFunc::func( Real const z ) const
{
	Real deriv( 0.0 );
	return func( z, false /*calc_deriv*/, deriv );
}

Real
GaussianChainGeneralFunc::dfunc( Real const z ) const
{
	Real deriv( 0.0 );
	func( z, true /*calc_deriv*/ , deriv );
	return deriv;
}

void
GaussianChainGeneralFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	Real dist;
	if ( in.good() ) {
		in >> dist;
		other_distances_.push_back( dist );
	}
	initialize_parameters();
}

void
GaussianChainGeneralFunc::show_definition( std::ostream &out ) const {
	out << "GAUSS_CHAIN " <<
		' ' << loop_fixed_cost_ <<
		' ' << gaussian_variance_;
	for ( auto distance : other_distances_ ) out << distance;
	out << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainGeneralFunc::GaussianChainGeneralFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainGeneralFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( other_distances_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_total_ ) ); // Real
	arc( CEREAL_NVP( N_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainGeneralFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( other_distances_ ); // utility::vector1<Real>
	arc( kB_T_ ); // Real
	arc( loop_fixed_cost_total_ ); // Real
	arc( N_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainGeneralFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainGeneralFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainGeneralFunc )
#endif // SERIALIZATION
