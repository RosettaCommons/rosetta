// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/ScalarWeightedFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


// C++ Headers

#include <sstream>
#include <iostream>
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

ScalarWeightedFunc::~ScalarWeightedFunc() = default;

bool ScalarWeightedFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< ScalarWeightedFunc const & > (other) );
	if ( weight_ != other_downcast.weight_ ) return false;

	return func_to_weight_ == other_downcast.func_to_weight_ ||
		( func_to_weight_ && other_downcast.func_to_weight_ &&
		*func_to_weight_ == *other_downcast.func_to_weight_ );
}

bool ScalarWeightedFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< ScalarWeightedFunc const * > ( &other );
}

Real
ScalarWeightedFunc::func( Real const x ) const
{
	return func_to_weight_->func( x ) * weight_;
}


Real
ScalarWeightedFunc::dfunc( Real const x ) const
{
	return func_to_weight_->dfunc( x ) * weight_;
}

void
ScalarWeightedFunc::read_data( std::istream& in )
{
	in >> weight_;

	FuncFactory func_factory;
	std::string func_type;
	in >> func_type;
	func_to_weight_ = func_factory.new_func( func_type );
	func_to_weight_->read_data( in );
}


void
ScalarWeightedFunc::show_definition( std::ostream &out ) const
{
	out << "SCALARWEIGHTEDFUNC  " << weight_ << " ";
	func_to_weight_->show_definition(out);
}

Size
ScalarWeightedFunc::show_violations(std::ostream &out, Real x, Size verbose_level, Real threshold) const
{
	out << "SCALARWEIGHTEDFUNC with weight:  " << weight_ << std::endl;
	func_to_weight_->show_definition(out);
	out << " with verbose_level " << verbose_level << ", threshold " << threshold << " and weighted_score " << ScalarWeightedFunc::func(x) << std::endl;
	return func_to_weight_->show_violations(out,x,verbose_level,threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::ScalarWeightedFunc::ScalarWeightedFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::ScalarWeightedFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( func_to_weight_ ) ); // FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::ScalarWeightedFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( weight_ ); // Real
	arc( func_to_weight_ ); // FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::ScalarWeightedFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::ScalarWeightedFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_ScalarWeightedFunc )
#endif // SERIALIZATION
