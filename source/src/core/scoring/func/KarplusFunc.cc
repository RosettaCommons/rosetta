// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/KarplusFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author Nikolas Sgourakis

#include <core/scoring/func/KarplusFunc.hh>
#include <core/types.hh>

#include <iostream>

#include <cmath>


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

bool KarplusFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< KarplusFunc const & > (other) );
	if ( A_      != other_downcast.A_      ) return false;
	if ( B_      != other_downcast.B_      ) return false;
	if ( C_      != other_downcast.C_      ) return false;
	if ( Dphi_   != other_downcast.Dphi_   ) return false;
	if ( x0_     != other_downcast.x0_     ) return false;
	if ( sd_     != other_downcast.sd_     ) return false;
	if ( offset_ != other_downcast.offset_ ) return false;

	return true;
}

bool KarplusFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< KarplusFunc const * > ( &other );
}

Real
KarplusFunc::func( Real const x ) const {
	Real const cosine = cos( x + Dphi_ );
	Real const j = A_ * cosine * cosine + B_ * cosine + C_;
	Real const z = (j - x0_)  / sd_;
	return z * z;

}

Real
KarplusFunc::dfunc( Real const x ) const {
	Real const sine = sin (x+Dphi_);
	Real const cosine = cos (x+Dphi_);
	Real const j = A_ * cosine * cosine + B_ * cosine + C_;
	return  -2 *( (j -x0_ ) / (sd_ * sd_) ) * sine *(2*A_*cosine + B_ ) ;
}

void
KarplusFunc::read_data( std::istream & in ) {
	in >> A_ >> B_ >>C_ >> Dphi_ >> x0_ >> sd_ >>offset_;
}

void KarplusFunc::show_definition( std::ostream & out ) const
{
	out << "KarplusFunc " << A_ <<' ' <<  B_ << ' ' << C_ << ' ' << Dphi_ << ' '<< x0_ << ' ' << sd_ << ' ' << offset_;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::KarplusFunc::KarplusFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::KarplusFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( A_ ) ); // Real
	arc( CEREAL_NVP( B_ ) ); // Real
	arc( CEREAL_NVP( C_ ) ); // Real
	arc( CEREAL_NVP( Dphi_ ) ); // Real
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( sd_ ) ); // Real
	arc( CEREAL_NVP( offset_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::KarplusFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( A_ ); // Real
	arc( B_ ); // Real
	arc( C_ ); // Real
	arc( Dphi_ ); // Real
	arc( x0_ ); // Real
	arc( sd_ ); // Real
	arc( offset_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::KarplusFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::KarplusFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_KarplusFunc )
#endif // SERIALIZATION
