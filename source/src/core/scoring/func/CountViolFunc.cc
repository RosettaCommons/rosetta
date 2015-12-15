// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CountViolFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/CountViolFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>

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

bool CountViolFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	CountViolFunc const & other_downcast( static_cast< CountViolFunc const & > (other) );
	if ( weight_      != other_downcast.weight_      ) return false;

	return func_to_weight_ == other_downcast.func_to_weight_ ||
		( func_to_weight_ && other_downcast.func_to_weight_ &&
		*func_to_weight_ == *other_downcast.func_to_weight_ );

	return true;
}

bool CountViolFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< CountViolFunc const * > ( &other );
}

Real
CountViolFunc::func( Real const x ) const
{
	return func_to_weight_->func( x ) * weight_;
}

Real
CountViolFunc::dfunc( Real const x ) const
{
	return func_to_weight_->dfunc( x ) * weight_;
}

void
CountViolFunc::read_data( std::istream & in )
{
	in >> count_viols_;
	weight_ = 1.0;
	FuncFactory func_factory;
	std::string func_type;
	in >> func_type;
	func_to_weight_ = func_factory.new_func( func_type );
	func_to_weight_->read_data( in );
}

/// @brief show some sort of stringified representation of the violations for this constraint.
core::Size CountViolFunc::show_violations(
	std::ostream& out,
	Real r,
	Size verbose_level,
	Real threshold
) const
{
	Size ct ( func_to_weight_->show_violations( out, r, verbose_level, threshold ) );
	// if ( verbose_level == 100 ) {
	count_viols_ += ct;
	// }
	return ct;
}

void
CountViolFunc::show_definition( std::ostream &out ) const
{
	using namespace ObjexxFCL::format;
	out << "COUNTVIOLFUNC" << RJ( 7, count_viols_ ) << " ";
	func_to_weight_->show_definition( out );
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::CountViolFunc::CountViolFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::CountViolFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( count_viols_ ) ); // core::Size
	arc( CEREAL_NVP( func_to_weight_ ) ); // FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::CountViolFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( weight_ ); // Real
	arc( count_viols_ ); // core::Size
	arc( func_to_weight_ ); // FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::CountViolFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::CountViolFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_CountViolFunc )
#endif // SERIALIZATION
