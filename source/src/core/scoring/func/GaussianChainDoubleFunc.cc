// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainDoubleFunc.cc
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainDoubleFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <numeric/constants.hh>

// C++ Headers
#include <iostream>
#include <cmath>

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

GaussianChainDoubleFunc::GaussianChainDoubleFunc( Real const gaussian_variance,
	Real const loop_fixed_cost,
	Real const D2 ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost ),
	D2_( D2 )
{
	initialize_parameters();
}

FuncOP
GaussianChainDoubleFunc::clone() const
{
	return FuncOP( new GaussianChainDoubleFunc( gaussian_variance_, loop_fixed_cost_, D2_ ) );
}

bool GaussianChainDoubleFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< GaussianChainDoubleFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( D2_                    != other_downcast.D2_                    ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( loop_fixed_cost_total_ != other_downcast.loop_fixed_cost_total_ ) return false;

	return true;
}

bool GaussianChainDoubleFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainDoubleFunc const * > ( &other );
}

void
GaussianChainDoubleFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	recompute_parameters();
}

void
GaussianChainDoubleFunc::recompute_parameters(){
	// this is a 'prefactor' in the probability.
	loop_fixed_cost_total_ = loop_fixed_cost_ + 1.5 * kB_T_ * log( Real(2.0 * pi) );
	loop_fixed_cost_total_ += kB_T_ * log( 2.0 );
	loop_fixed_cost_total_ += 0.5 * kB_T_ * log( gaussian_variance_ );
	// further normalization factor is put into func itself.
}

Real
GaussianChainDoubleFunc::func( Real const z ) const
{
	Real const & D1 = z;
	Real const & D2 = D2_;

	Real const term1 = exp( -( D1 - D2) * (D1 - D2)/ (2 * gaussian_variance_ ) );
	Real const term2 = exp( -( D1 + D2) * (D1 + D2)/ (2 * gaussian_variance_ ) );
	Real const loop_energy = -kB_T_ * log( ( term1 - term2 ) / (D1 * D2 ));
	return ( loop_fixed_cost_total_ + loop_energy );
}

Real
GaussianChainDoubleFunc::dfunc( Real const z ) const
{
	Real const & D1 = z;
	Real const & D2 = D2_;

	/////
	Real const first_logderiv_term = kB_T_ / z;

	/////
	Real const term1 = exp( -( D1 - D2) * (D1 - D2)/ (2 * gaussian_variance_ ) );
	Real const term2 = exp( -( D1 + D2) * (D1 + D2)/ (2 * gaussian_variance_ ) );

	Real second_logderiv_term = ( ( D2 - D1 ) * term1 + (D1 + D2) * term2 )/ (term1 - term2 );
	second_logderiv_term *= ( -kB_T_ / gaussian_variance_ );

	return first_logderiv_term + second_logderiv_term;
}

void
GaussianChainDoubleFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	initialize_parameters();
}

void
GaussianChainDoubleFunc::show_definition( std::ostream &out ) const {
	out << "GAUSS_CHAIN " <<
		' ' << loop_fixed_cost_ <<
		' ' << gaussian_variance_ << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainDoubleFunc::GaussianChainDoubleFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainDoubleFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( D2_ ) ); // Real
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_total_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainDoubleFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( D2_ ); // Real
	arc( kB_T_ ); // Real
	arc( loop_fixed_cost_total_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainDoubleFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainDoubleFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainDoubleFunc )
#endif // SERIALIZATION
