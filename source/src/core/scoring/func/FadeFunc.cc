// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/FadeFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das


#include <core/scoring/func/FadeFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>


// C++ Headers
#include <iostream>
#include <sstream>

#include <utility/vector1.hh>


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

FuncOP
FadeFunc::clone() const
{
	return FuncOP( new FadeFunc( cutoff_lower_, cutoff_upper_,
		fade_zone_, well_depth_,
		well_offset_ ) );
}

bool FadeFunc::operator == ( Func const & other ) const
{
	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< FadeFunc const & > (other) );
	if ( cutoff_lower_ != other_downcast.cutoff_lower_ ) return false;
	if ( cutoff_upper_ != other_downcast.cutoff_upper_ ) return false;
	if ( fade_zone_    != other_downcast.fade_zone_    ) return false;
	if ( well_depth_   != other_downcast.well_depth_   ) return false;
	if ( well_offset_  != other_downcast.well_offset_  ) return false;
	return true;
}

bool FadeFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< FadeFunc const * > ( &other );
}

Real
FadeFunc::func( Real const z ) const
{
	Real fade_value( 1.0 );

	if ( z < cutoff_lower_ || z > cutoff_upper_ ) {
		fade_value = 0.0;
	} else if ( z < cutoff_lower_ + fade_zone_ ) {
		//Check little strip near lower cutoff.
		Real const b = -1.0 * ( z - (cutoff_lower_ + fade_zone_) )/ fade_zone_;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		//  fade_deriv = -1.0 * (6 * b2 - 6 * b ) / fade_zone_;
	} else if ( z > cutoff_upper_ - fade_zone_ ) {
		//Check little strip near upper cutoff.
		Real const b =  ( z - (cutoff_upper_ - fade_zone_) )/ fade_zone_;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		//  fade_deriv = (6 * b2 - 6 * b ) / fade_zone_;
	}

	return well_depth_ * fade_value + well_offset_;
}

Real
FadeFunc::dfunc( Real const z ) const
{
	Real fade_deriv( 0.0 );

	if ( z < cutoff_lower_ || z > cutoff_upper_ ) {
		fade_deriv = 0.0;
	} else if ( z < cutoff_lower_ + fade_zone_ ) {
		//Check little strip near lower cutoff.
		Real const b = -1.0 * ( z - (cutoff_lower_ + fade_zone_) )/ fade_zone_;
		Real const b2 = b*b;
		fade_deriv = -1.0 * (6 * b2 - 6 * b ) / fade_zone_;
	} else if ( z > cutoff_upper_ - fade_zone_ ) {
		//Check little strip near upper cutoff.
		Real const b =  ( z - (cutoff_upper_ - fade_zone_) )/ fade_zone_;
		Real const b2 = b*b;
		fade_deriv = (6 * b2 - 6 * b ) / fade_zone_;
	}

	return well_depth_ * fade_deriv;
}

void
FadeFunc::read_data( std::istream& in ) {

	in >> cutoff_lower_;
	in >> cutoff_upper_;
	in >> fade_zone_;
	in >> well_depth_ ;

	well_offset_ = 0.0;

	// there may be another number, which would be the well offset.
	char dummy;
	while ( in.good() && ( in.peek() == ' ' ) ) in.get( dummy );
	if ( in.good() ) {
		char c = in.peek();
		if ( isdigit( c ) || c == '-' || c == '+' ) { // allow for negative well offsets...
			in >> well_offset_;
		}
	}
}

void
FadeFunc::show_definition( std::ostream &out ) const {
	out << "FADE " <<
		' ' << cutoff_lower_ <<
		' ' << cutoff_upper_ <<
		' ' << fade_zone_ <<
		' ' << well_depth_ <<
		' ' << well_offset_ << std::endl;
}

Size
FadeFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "FADE " <<  ( x < cutoff_lower_ || x > cutoff_upper_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::FadeFunc::FadeFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::FadeFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( cutoff_lower_ ) ); // Real
	arc( CEREAL_NVP( cutoff_upper_ ) ); // Real
	arc( CEREAL_NVP( fade_zone_ ) ); // Real
	arc( CEREAL_NVP( well_depth_ ) ); // Real
	arc( CEREAL_NVP( well_offset_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::FadeFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( cutoff_lower_ ); // Real
	arc( cutoff_upper_ ); // Real
	arc( fade_zone_ ); // Real
	arc( well_depth_ ); // Real
	arc( well_offset_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::FadeFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::FadeFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_FadeFunc )
#endif // SERIALIZATION
