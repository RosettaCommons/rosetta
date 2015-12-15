// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularSigmoidalFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author Robert Vernon

#include <core/scoring/func/CircularSigmoidalFunc.hh>
#include <numeric/angle.functions.hh>
#include <core/types.hh>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif
#include <cmath>

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

bool CircularSigmoidalFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	CircularSigmoidalFunc const & other_downcast( static_cast< CircularSigmoidalFunc const & > (other) );
	if ( xC_     != other_downcast.xC_     ) return false;
	if ( m_      != other_downcast.m_      ) return false;
	if ( o1_     != other_downcast.o1_     ) return false;
	if ( o2_     != other_downcast.o2_     ) return false;
	if ( offset_ != other_downcast.offset_ ) return false;
	return true;
}

bool CircularSigmoidalFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< CircularSigmoidalFunc const * > ( &other );
}

Real
CircularSigmoidalFunc::func( Real const x ) const {
	Real const x0 = numeric::nearest_angle_radians(x,xC_)-xC_;
	Real const z = offset_ + (1/(1+ std::pow( M_E, (-m_*(x0-o1_)) ))) - (1/(1+ std::pow( M_E, (-m_*(x0-o2_)) )));

	return z;
}

Real
CircularSigmoidalFunc::dfunc( Real const x ) const {

	Real const x0 = numeric::nearest_angle_radians(x,xC_)-xC_;

	Real const z = m_*std::pow(M_E,m_*x0+m_*o1_)/(2*std::pow(M_E,m_*x0+m_*o1_)+std::pow(M_E,2*m_*x0)+std::pow(M_E,2*m_*o1_))
		- m_*std::pow(M_E,m_*x0+m_*o2_)/(2*std::pow(M_E,m_*x0+m_*o2_)+std::pow(M_E,2*m_*x0)+std::pow(M_E,2*m_*o2_));

	return z;
}

void
CircularSigmoidalFunc::read_data( std::istream & in ) {
	in >> xC_ >> m_ >> o1_ >> o2_;
}

void CircularSigmoidalFunc::show_definition( std::ostream & out ) const
{
	out << "CircularSigmoidalFunc " << xC_ << ' ' << m_ << ' ' << o1_ << ' ' << o2_;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::CircularSigmoidalFunc::CircularSigmoidalFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::CircularSigmoidalFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( xC_ ) ); // Real
	arc( CEREAL_NVP( m_ ) ); // Real
	arc( CEREAL_NVP( o1_ ) ); // Real
	arc( CEREAL_NVP( o2_ ) ); // Real
	arc( CEREAL_NVP( offset_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::CircularSigmoidalFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( xC_ ); // Real
	arc( m_ ); // Real
	arc( o1_ ); // Real
	arc( o2_ ); // Real
	arc( offset_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::CircularSigmoidalFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::CircularSigmoidalFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_CircularSigmoidalFunc )
#endif // SERIALIZATION
