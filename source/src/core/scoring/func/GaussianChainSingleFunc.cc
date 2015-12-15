// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/GaussianChainSingleFunc.hh
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainSingleFunc.hh>
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

GaussianChainSingleFunc::GaussianChainSingleFunc( Real const gaussian_variance,
	Real const loop_fixed_cost ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost )
{
	initialize_parameters();
}

FuncOP
GaussianChainSingleFunc::clone() const
{
	return FuncOP( new GaussianChainSingleFunc( gaussian_variance_, loop_fixed_cost_ ) );
}

bool GaussianChainSingleFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	GaussianChainSingleFunc const & other_downcast( static_cast< GaussianChainSingleFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( loop_fixed_cost_total_ != other_downcast.loop_fixed_cost_total_ ) return false;
	return true;
}

bool GaussianChainSingleFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainSingleFunc const * > ( &other );
}


void
GaussianChainSingleFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	recompute_parameters();
}

void
GaussianChainSingleFunc::recompute_parameters(){
	// the additional term is the normalization factor for a Gaussian distribution.
	loop_fixed_cost_total_ = loop_fixed_cost_ + 1.5 * kB_T_ * log( (2 * pi * gaussian_variance_ ) );
}

Real
GaussianChainSingleFunc::func( Real const z ) const
{
	Real const loop_harmonic = kB_T_ * (z * z)/ ( 2 * gaussian_variance_ );
	return ( loop_fixed_cost_total_ + loop_harmonic );
}

Real
GaussianChainSingleFunc::dfunc( Real const z ) const
{
	Real const loop_harmonic_deriv = kB_T_ * z / gaussian_variance_;
	return loop_harmonic_deriv;
}

void
GaussianChainSingleFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	initialize_parameters();
}

void
GaussianChainSingleFunc::show_definition( std::ostream &out ) const {
	out << "GAUSS_CHAIN " <<
		' ' << loop_fixed_cost_ <<
		' ' << gaussian_variance_ << std::endl;
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainSingleFunc::GaussianChainSingleFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainSingleFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_total_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainSingleFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( kB_T_ ); // Real
	arc( loop_fixed_cost_total_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainSingleFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainSingleFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainSingleFunc )
#endif // SERIALIZATION
