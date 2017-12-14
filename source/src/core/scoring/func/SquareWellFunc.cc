// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SquareWellFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das


#include <core/scoring/func/SquareWellFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


#include <utility/vector1.hh>
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

bool SquareWellFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< SquareWellFunc const & > (other) );
	if ( x0_         != other_downcast.x0_         ) return false;
	if ( well_depth_ != other_downcast.well_depth_ ) return false;
	return true;
}

bool SquareWellFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SquareWellFunc const * > ( &other );
}

Real
SquareWellFunc::func( Real const x ) const
{
	if ( x < x0_ ) {
		return well_depth_;
	}
	return 0.0;
}

Real
SquareWellFunc::dfunc( Real const /*x*/ ) const
{
	return 0.0; //This is bad news for the minimizer...
}

void
SquareWellFunc::read_data( std::istream& in ) {
	in >> x0_ >> well_depth_;
}

void
SquareWellFunc::show_definition( std::ostream &out ) const {
	out << "SQUARE_WELL " << x0_ << " " << well_depth_ << std::endl;
}

Size
SquareWellFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "SQUARE_WELL " <<  ( x < x0_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::SquareWellFunc::SquareWellFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SquareWellFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( well_depth_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SquareWellFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( well_depth_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SquareWellFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::SquareWellFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SquareWellFunc )
#endif // SERIALIZATION
