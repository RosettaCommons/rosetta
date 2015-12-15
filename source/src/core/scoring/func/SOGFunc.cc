// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SOGFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/SOGFunc_Impl.hh>
#include <core/types.hh>

// C++ Headers

#include <utility/vector1.hh>


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

SOGFunc::SOGFunc(
	const utility::vector1< core::Real >& means,
	const utility::vector1< core::Real >& sdevs,
	const utility::vector1< core::Real >& weights
) : member_(means, sdevs, weights) {}

SOGFunc::SOGFunc(
	core::Real mean,
	core::Real sdev
) : member_(utility::vector1< core::Real >(1,mean),
	utility::vector1< core::Real >(1,sdev),
	utility::vector1< core::Real >(1,1.0)) {}


bool SOGFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	SOGFunc const & other_downcast( static_cast< SOGFunc const & > (other) );
	return member_ == other_downcast.member_;

}

bool SOGFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SOGFunc const * > ( &other );
}

void
SOGFunc::read_data( std::istream & in ) {
	member_.read_data( in );
}

void
SOGFunc::clear_() {
	member_.clear_();
}

core::Real
SOGFunc::get_alt_score_( Real const x ) const {
	return member_.get_alt_score_(x);
}

Real
SOGFunc::func( Real const x ) const {
	return member_.func(x);
} // func

Real
SOGFunc::dfunc( Real const x ) const {
	return member_.dfunc(x);
} // dfunc

void SOGFunc::check_bounds( Real const x, Real const val ) const {
	member_.check_bounds(x,val);
}

void SOGFunc::show_definition( std::ostream & out ) const {
	member_.show_definition(out);
} // show_definition

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SOGFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( member_ ) ); // class core::scoring::func::SOGFunc_Impl
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SOGFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( member_ ); // class core::scoring::func::SOGFunc_Impl
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SOGFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::SOGFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SOGFunc )
#endif // SERIALIZATION
