// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/ConstantFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author John Karanicolas

#include <core/scoring/func/ConstantFunc.hh>
#include <core/types.hh>

// C++ Headers
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

bool ConstantFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< ConstantFunc const & > (other) );
	return return_val_ == other_downcast.return_val_;
}

bool ConstantFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< ConstantFunc const * > ( &other );
}

void
ConstantFunc::read_data( std::istream & in ) {
	in  >> return_val_;
}

Real
ConstantFunc::func( Real const ) const {
	return return_val_;
} // func

Real
ConstantFunc::dfunc( Real const ) const {
	return 0;
} // dfunc

void ConstantFunc::show_definition( std::ostream & out ) const {
	out << "CONSTANTFUNC " << return_val_ << "\n";
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::ConstantFunc::ConstantFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::ConstantFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( return_val_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::ConstantFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( return_val_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::ConstantFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::ConstantFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_ConstantFunc )
#endif // SERIALIZATION
