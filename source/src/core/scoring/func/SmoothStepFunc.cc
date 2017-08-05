// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SmoothStepFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Joaquin Ambia, Jason K. Lai


#include <core/scoring/func/SmoothStepFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <cmath>

#include <utility/vector1.hh>
#include <sstream>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION



// C++ Headers


namespace core {
namespace scoring {
namespace func {

bool SmoothStepFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	SmoothStepFunc const & other_downcast( static_cast< SmoothStepFunc const & > (other) );
	if ( low_         != other_downcast.low_         ) return false;
	if ( high_        != other_downcast.high_        ) return false;
	return true;
}

bool SmoothStepFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SmoothStepFunc const * > ( &other );
}

void
SmoothStepFunc::read_data( std::istream & in ) {
	in >> low_ >> high_;
}


Real
SmoothStepFunc::func( Real const x ) const
{
	if ( x < low_ ) {
		return 0;
	} else if ( x > high_ ) {
		return 1;
	} else {
		Real a = (x - low_)/(high_ - low_);
		return 6*std::pow(a, 5) - 15*std::pow(a, 4) + 10*std::pow(a,3);
	}
}

Real
SmoothStepFunc::dfunc( Real const x ) const
{
	if ( x < low_ || x > high_ ) {
		return 0;
	} else {
		Real a = (x - low_)/(high_ - low_);
		return (30*std::pow(a, 4) - 60*std::pow(a, 3) + 30*std::pow(a,2))/(high_ - low_);
	}
}

void
SmoothStepFunc::show_definition( std::ostream &out ) const {
	out << low_ << high_ << std::endl;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::SmoothStepFunc::SmoothStepFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SmoothStepFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( low_ ) ); // Real
	arc( CEREAL_NVP( high_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SmoothStepFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( low_ ); // Real
	arc( high_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SmoothStepFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::SmoothStepFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SmoothStepFunc )
#endif // SERIALIZATION
