// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/LinearPenaltyFunction.hh
/// @author Dominik Gront


#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/types.hh>

// C++ Headers
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

bool LinearPenaltyFunction::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	LinearPenaltyFunction const & other_downcast( static_cast< LinearPenaltyFunction const & > (other) );
	if ( x_middle_   != other_downcast.x_middle_   ) return false;
	if ( well_depth_ != other_downcast.well_depth_ ) return false;
	if ( half_width_ != other_downcast.half_width_ ) return false;
	if ( slope_      != other_downcast.slope_      ) return false;

	return true;
}

bool LinearPenaltyFunction::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< LinearPenaltyFunction const * > ( &other );
}

Real
LinearPenaltyFunction::func( Real const x ) const {

	Real dev = fabs(x-x_middle_);
	if ( dev <= half_width_ ) {
		return well_depth_;
	}
	return well_depth_ + slope_ * (dev-half_width_);
}

Real
LinearPenaltyFunction::dfunc( Real const /*x*/ ) const
{
	return 0.0;
}

void
LinearPenaltyFunction::read_data( std::istream& in ) {
	in >> x_middle_ >> well_depth_ >> half_width_ >> slope_;
}

void
LinearPenaltyFunction::show_definition( std::ostream &out ) const {
	out << "LINEAR_PENALTY " << x_middle_ << " " << well_depth_ <<  " " << half_width_ <<
		" " << slope_ << std::endl;
}

Size
LinearPenaltyFunction::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {

	if ( verbose_level > 100 ) {
		out << "LINEAR_PENALTY " <<  ( x < x_middle_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::LinearPenaltyFunction::LinearPenaltyFunction() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::LinearPenaltyFunction::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x_middle_ ) ); // Real
	arc( CEREAL_NVP( well_depth_ ) ); // Real
	arc( CEREAL_NVP( half_width_ ) ); // Real
	arc( CEREAL_NVP( slope_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::LinearPenaltyFunction::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x_middle_ ); // Real
	arc( well_depth_ ); // Real
	arc( half_width_ ); // Real
	arc( slope_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::LinearPenaltyFunction );
CEREAL_REGISTER_TYPE( core::scoring::func::LinearPenaltyFunction )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_LinearPenaltyFunction )
#endif // SERIALIZATION
