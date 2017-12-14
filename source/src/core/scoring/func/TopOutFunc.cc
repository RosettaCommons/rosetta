// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/TopOutFunc.hh
/// @brief Implementation of phenix "top-out" function
///   Similar to Geman-McClure: harmonic near 'x0_', flat past 'limit_'
/// @author Frank DiMaio


#include <core/scoring/func/TopOutFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <cmath>
#include <sstream>

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

TopOutFunc::TopOutFunc( Real weight_in, Real x0_in, Real limit_in ) :
	x0_( x0_in ), weight_( weight_in), limit_( limit_in ) {}

FuncOP TopOutFunc::clone() const { return FuncOP( new TopOutFunc( *this ) ); }

bool TopOutFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< TopOutFunc const & > (other) );
	if ( x0_ != other_downcast.x0_ ) return false;
	if ( weight_ != other_downcast.weight_ ) return false;
	if ( limit_ != other_downcast.limit_ ) return false;

	return true;
}

bool TopOutFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< TopOutFunc const * > ( &other );
}

Real
TopOutFunc::func( Real const x ) const {
	Real xoff = x-x0_;
	Real top = weight_ * limit_ * limit_;
	Real score = top * (1.0 - exp(-weight_*xoff*xoff/top));
	return score;
}

Real
TopOutFunc::dfunc( Real const x ) const {
	Real xoff = x-x0_;
	Real top = weight_ * limit_ * limit_;
	Real grad = 2.0 * weight_ * xoff * exp(-(weight_*xoff*xoff)/top);
	return (grad);
}

void
TopOutFunc::read_data( std::istream& in ) {
	in >> weight_ >> x0_ >> limit_;
}

void
TopOutFunc::show_definition( std::ostream &out ) const {
	out << "TOPOUT " << weight_ << " " << x0_ << " " << limit_ << std::endl;
}

Size
TopOutFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "TOPOUT " <<  func(x) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::TopOutFunc::TopOutFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::TopOutFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( limit_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::TopOutFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( weight_ ); // Real
	arc( limit_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::TopOutFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::TopOutFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_TopOutFunc )
#endif // SERIALIZATION
