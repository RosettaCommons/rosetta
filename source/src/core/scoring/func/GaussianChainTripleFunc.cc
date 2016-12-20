// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainTripleFunc.hh
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainTripleFunc.hh>
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

GaussianChainTripleFunc::GaussianChainTripleFunc( Real const gaussian_variance,
	Real const loop_fixed_cost,
	Real const D2, Real const D3 ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost ),
	D2_( D2 ),
	D3_( D3 )
{
	initialize_parameters();
}

FuncOP
GaussianChainTripleFunc::clone() const
{
	return FuncOP( new GaussianChainTripleFunc( gaussian_variance_, loop_fixed_cost_, D2_, D3_ ) );
}

bool GaussianChainTripleFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	GaussianChainTripleFunc const & other_downcast( static_cast< GaussianChainTripleFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( D2_                    != other_downcast.D2_                    ) return false;
	if ( D3_                    != other_downcast.D3_                    ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( loop_fixed_cost_total_ != other_downcast.loop_fixed_cost_total_ ) return false;

	return true;
}

bool GaussianChainTripleFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainTripleFunc const * > ( &other );
}

void
GaussianChainTripleFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	recompute_parameters();
}

void
GaussianChainTripleFunc::recompute_parameters(){
	// this is a 'prefactor' in the probability.
	loop_fixed_cost_total_ = loop_fixed_cost_ + kB_T_ * log( 16 * pi );
	// further normalization factor is put into func itself.
}

Real
GaussianChainTripleFunc::func( Real const z ) const
{
	Real const & D1 = z;
	Real const & D2 = D2_;
	Real const & D3 = D3_;

	Real const s = sqrt( 2 * gaussian_variance_ );

	// yea, error functions! How crazy is that.
#ifdef WIN32
	Real const term0 = boost::math::erf( ( D1 + D2 + D3 )/ s );
	Real const term1 = boost::math::erf( (-D1 + D2 + D3 )/ s );
	Real const term2 = boost::math::erf( ( D1 - D2 + D3 )/ s );
	Real const term3 = boost::math::erf( ( D1 + D2 - D3 )/ s );
#else
	Real const term0 = erf( ( D1 + D2 + D3 )/ s );
	Real const term1 = erf( (-D1 + D2 + D3 )/ s );
	Real const term2 = erf( ( D1 - D2 + D3 )/ s );
	Real const term3 = erf( ( D1 + D2 - D3 )/ s );
#endif
	Real const loop_energy = -kB_T_ * log( ( -term0 + term1 + term2 + term3 ) / (D1 * D2 * D3));

	return ( loop_fixed_cost_total_ + loop_energy );
}

Real
GaussianChainTripleFunc::dfunc( Real const z ) const
{

	Real const & D1 = z;
	Real const & D2 = D2_;
	Real const & D3 = D3_;

	/////
	Real const first_logderiv_term = kB_T_ / z;

	Real const s = sqrt( 2 * gaussian_variance_ );

#ifdef WIN32
	Real const term0 = boost::math::erf( ( D1 + D2 + D3 )/ s );
	Real const term1 = boost::math::erf( (-D1 + D2 + D3 )/ s );
	Real const term2 = boost::math::erf( ( D1 - D2 + D3 )/ s );
	Real const term3 = boost::math::erf( ( D1 + D2 - D3 )/ s );
#else
	Real const term0 = erf( ( D1 + D2 + D3 )/ s );
	Real const term1 = erf( (-D1 + D2 + D3 )/ s );
	Real const term2 = erf( ( D1 - D2 + D3 )/ s );
	Real const term3 = erf( ( D1 + D2 - D3 )/ s );
#endif
	Real const deriv_term0 = exp( -1.0 * pow( ( D1 + D2 + D3 )/ s, 2 ) );
	Real const deriv_term1 = exp( -1.0 * pow( (-D1 + D2 + D3 )/ s, 2 ) );
	Real const deriv_term2 = exp( -1.0 * pow( ( D1 - D2 + D3 )/ s, 2 ) );
	Real const deriv_term3 = exp( -1.0 * pow( ( D1 + D2 - D3 )/ s, 2 ) );

	Real second_logderiv_term = -deriv_term0  - deriv_term1 + deriv_term2 + deriv_term3;
	second_logderiv_term /= ( -term0 + term1 + term2 + term3 );
	second_logderiv_term *= -1.0 * kB_T_ * 2 / ( sqrt( pi ) * s );

	return first_logderiv_term + second_logderiv_term;

}

void
GaussianChainTripleFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	initialize_parameters();
}

void
GaussianChainTripleFunc::show_definition( std::ostream &out ) const {
	out << "GAUSS_CHAIN " <<
		' ' << loop_fixed_cost_ <<
		' ' << gaussian_variance_ << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainTripleFunc::GaussianChainTripleFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainTripleFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( D2_ ) ); // Real
	arc( CEREAL_NVP( D3_ ) ); // Real
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_total_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainTripleFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( D2_ ); // Real
	arc( D3_ ); // Real
	arc( kB_T_ ); // Real
	arc( loop_fixed_cost_total_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainTripleFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainTripleFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainTripleFunc )
#endif // SERIALIZATION
