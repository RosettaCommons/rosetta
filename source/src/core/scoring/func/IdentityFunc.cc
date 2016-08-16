// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/IdentityFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/func/IdentityFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <sstream>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

bool IdentityFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	return true;
}

bool IdentityFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< IdentityFunc const * > ( &other );
}

Real
IdentityFunc::func( Real const x ) const {
	return x;
}

Real
IdentityFunc::dfunc( Real const /*x*/ ) const {
	return 1;
}

void
IdentityFunc::read_data( std::istream & /*in*/ ) {}

void
IdentityFunc::show_definition( std::ostream & out ) const {
	out << "IDENTITY";
}

Size
IdentityFunc::show_violations(
	std::ostream & out,
	Real x,
	Size verbose_level,
	Real threshold
) const {
	return Func::show_violations( out, x, verbose_level, threshold );
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::IdentityFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::IdentityFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::IdentityFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::IdentityFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_IdentityFunc )
#endif // SERIALIZATION
