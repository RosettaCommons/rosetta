// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/CircularPowerFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <core/scoring/func/CircularPowerFunc.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>
#include <utility/assert.hh>

// C++ Headers

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

CircularPowerFunc::CircularPowerFunc(
	Real const x0_radians,
	Real const sd_radians,
	int const power,
	Real const weight
) :
	x0_( x0_radians ),
	sd_( sd_radians ),
	power_( power ),
	weight_( weight )
{
	debug_assert( std::abs( std::pow(3.0,2) - 9.0 ) < 1e-3 );
	debug_assert( power_ != 1 && power_ != 0 );
}

FuncOP
CircularPowerFunc::clone() const {
	return FuncOP( new CircularPowerFunc( *this ) );
}

bool CircularPowerFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	CircularPowerFunc const & other_downcast( static_cast< CircularPowerFunc const & > (other) );
	if ( x0_     != other_downcast.x0_     ) return false;
	if ( sd_     != other_downcast.sd_     ) return false;
	if ( power_  != other_downcast.power_  ) return false;
	if ( weight_ != other_downcast.weight_ ) return false;

	return true;
}

bool CircularPowerFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< CircularPowerFunc const * > ( &other );
}

Real
CircularPowerFunc::func( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return weight_ * std::pow( z, power_ );
}

Real
CircularPowerFunc::dfunc( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return weight_ * power_ * std::pow( z, power_ - 1 ) / sd_;
}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::CircularPowerFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // const Real
	arc( CEREAL_NVP( sd_ ) ); // const Real
	arc( CEREAL_NVP( power_ ) ); // const int
	arc( CEREAL_NVP( weight_ ) ); // const Real
}

/// @brief Hand crafted deserialization method
template< class Archive >
void
core::scoring::func::CircularPowerFunc::load_and_construct( Archive & arc, cereal::construct< core::scoring::func::CircularPowerFunc > & construct ) {
	Real x0, sd, weight; int power;
	arc( x0, sd, power, weight );
	construct( x0, sd, power, weight );
	// EXEMPT x0_ sd_ power_ weight_
}

SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( core::scoring::func::CircularPowerFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::CircularPowerFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_CircularPowerFunc )
#endif // SERIALIZATION
