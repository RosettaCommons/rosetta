// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainQuadrupleFunc.hh
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainQuadrupleFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <numeric/constants.hh>

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

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

GaussianChainQuadrupleFunc::GaussianChainQuadrupleFunc( Real const gaussian_variance,
	Real const loop_fixed_cost,
	Real const D2, Real const D3, Real const D4 ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost ),
	D2_( D2 ),
	D3_( D3 ),
	D4_( D4 )
{
	initialize_parameters();
}

FuncOP
GaussianChainQuadrupleFunc::clone() const
{
	return FuncOP( new GaussianChainQuadrupleFunc( gaussian_variance_, loop_fixed_cost_, D2_, D3_, D4_ ) );
}

bool GaussianChainQuadrupleFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	GaussianChainQuadrupleFunc const & other_downcast( static_cast< GaussianChainQuadrupleFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( D2_                    != other_downcast.D2_                    ) return false;
	if ( D3_                    != other_downcast.D3_                    ) return false;
	if ( D4_                    != other_downcast.D4_                    ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( loop_fixed_cost_total_ != other_downcast.loop_fixed_cost_total_ ) return false;

	return true;
}

bool GaussianChainQuadrupleFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainQuadrupleFunc const * > ( &other );
}

void
GaussianChainQuadrupleFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	recompute_parameters();
}

void
GaussianChainQuadrupleFunc::recompute_parameters(){
	// this is a 'prefactor' in the probability.
	loop_fixed_cost_total_ = loop_fixed_cost_;
	loop_fixed_cost_total_ += kB_T_ * log( 32.0 );
	loop_fixed_cost_total_ += 1.5 * kB_T_ * log( pi );
	loop_fixed_cost_total_ += -0.5 * kB_T_ * log( 2 * gaussian_variance_ );
	// further normalization factor is put into func itself.
}

///////////////////////////
Real
L( Real const & x ){
	// note that this is the integral of erf( x ).
	// up to a sqrt( pi ) factor.
	// that's the key to a general expression.
#ifdef WIN32
	return (sqrt( pi ) * x * boost::math::erf( x )) + exp( -1.0 * x * x);
#else
	return (sqrt( pi ) * x * erf( x )) + exp( -1.0 * x * x);
#endif
}

Real
GaussianChainQuadrupleFunc::func( Real const z ) const
{
	Real const & D1 = z;
	Real const & D2 = D2_;
	Real const & D3 = D3_;
	Real const & D4 = D4_;

	Real const s = sqrt( 2 * gaussian_variance_ );

	Real L_sum = 0.0;
	for ( int sgn2 = -1; sgn2 <= 1; sgn2 += 2 ) {
		for ( int sgn3 = -1; sgn3 <= 1; sgn3 += 2 ) {
			for ( int sgn4 = -1; sgn4 <= 1; sgn4 += 2 ) {
				Real D_sum = D1 + ( sgn2 * D2 ) + ( sgn3 * D3 ) + ( sgn4 * D4 );
				D_sum /= s;
				int sgn = -1 * sgn2 * sgn3 * sgn4;
				L_sum += sgn * L( D_sum );
			}
		}
	}

	Real const loop_energy = -kB_T_ * log( L_sum / (D1 * D2 * D3 * D4) );
	return ( loop_fixed_cost_total_ + loop_energy );
}

Real
GaussianChainQuadrupleFunc::dfunc( Real const z ) const
{
	Real const & D1 = z;
	Real const & D2 = D2_;
	Real const & D3 = D3_;
	Real const & D4 = D4_;

	/////
	Real const first_logderiv_term = kB_T_ / z;

	Real const s = sqrt( 2 * gaussian_variance_ );

	Real L_sum = 0.0;
	Real deriv_sum = 0.0;
	for ( int sgn2 = -1; sgn2 <= 1; sgn2 += 2 ) {
		for ( int sgn3 = -1; sgn3 <= 1; sgn3 += 2 ) {
			for ( int sgn4 = -1; sgn4 <= 1; sgn4 += 2 ) {
				Real D_sum = D1 + ( sgn2 * D2 ) + ( sgn3 * D3 ) + ( sgn4 * D4 );
				D_sum /= s;
				int sgn = -1 * sgn2 * sgn3 * sgn4;
				L_sum += sgn * L( D_sum );
#ifdef WIN32
				deriv_sum += sgn * sqrt( pi ) * boost::math::erf( D_sum );
#else
				deriv_sum += sgn * sqrt( pi ) * erf( D_sum );
#endif
			}
		}
	}

	Real second_logderiv_term = deriv_sum / L_sum;
	second_logderiv_term *= -1.0 * kB_T_ / s;

	return first_logderiv_term + second_logderiv_term;
}

void
GaussianChainQuadrupleFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	initialize_parameters();
}

void
GaussianChainQuadrupleFunc::show_definition( std::ostream &out ) const {
	out << "GAUSS_CHAIN " <<
		' ' << loop_fixed_cost_ <<
		' ' << gaussian_variance_ << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainQuadrupleFunc::GaussianChainQuadrupleFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainQuadrupleFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( D2_ ) ); // Real
	arc( CEREAL_NVP( D3_ ) ); // Real
	arc( CEREAL_NVP( D4_ ) ); // Real
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_total_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainQuadrupleFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( D2_ ); // Real
	arc( D3_ ); // Real
	arc( D4_ ); // Real
	arc( kB_T_ ); // Real
	arc( loop_fixed_cost_total_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainQuadrupleFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainQuadrupleFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainQuadrupleFunc )
#endif // SERIALIZATION
