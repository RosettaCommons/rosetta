// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/CircularHarmonicFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <numeric/angle.functions.hh>
#include <core/types.hh>

#include <iostream>

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

bool CircularHarmonicFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	CircularHarmonicFunc const & other_downcast( static_cast< CircularHarmonicFunc const & > (other) );
	if ( x0_     != other_downcast.x0_     ) return false;
	if ( sd_     != other_downcast.sd_     ) return false;
	if ( offset_ != other_downcast.offset_ ) return false;
	return true;
}

bool CircularHarmonicFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< CircularHarmonicFunc const * > ( &other );
}

Real
CircularHarmonicFunc::func( Real const x ) const {
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ )/sd_;
	return z * z + offset_;
}

Real
CircularHarmonicFunc::dfunc( Real const x ) const {
	return 2 * ( numeric::nearest_angle_radians(x,x0_)-x0_ ) / ( sd_ * sd_ );
}

void
CircularHarmonicFunc::read_data( std::istream & in ) {
	in >> x0_ >> sd_;
}

void CircularHarmonicFunc::show_definition( std::ostream & out ) const
{
	//out << "CircularHarmonicFunc " << x0_ << ' ' << sd_;
	out << "CIRCULARHARMONIC " << x0_ << ' ' << sd_ << std::endl;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::CircularHarmonicFunc::CircularHarmonicFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::CircularHarmonicFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( sd_ ) ); // Real
	arc( CEREAL_NVP( offset_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::CircularHarmonicFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( sd_ ); // Real
	arc( offset_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::CircularHarmonicFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::CircularHarmonicFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_CircularHarmonicFunc )
#endif // SERIALIZATION
