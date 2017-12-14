// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SquareWell2Func.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Rhiju Das
/// @author Jianqing Xu

#include <sstream>

#include <basic/Tracer.hh>
#include <core/scoring/func/SquareWell2Func.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

static basic::Tracer TR( "core.scoring.constraints.SquareWall2Func" );

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

bool SquareWell2Func::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< SquareWell2Func const & > (other) );
	if ( x0_         != other_downcast.x0_         ) return false;
	if ( x_range_    != other_downcast.x_range_    ) return false;
	if ( well_depth_ != other_downcast.well_depth_ ) return false;
	return true;
}

bool SquareWell2Func::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SquareWell2Func const * > ( &other );
}

Real
SquareWell2Func::func( Real const x ) const
{
	Real const z = ( numeric::nearest_angle_radians(x,x0_)-x0_ );
	// TR<<"x="<<x<<"   x0_="<<x0_<<"   |x-x0|="<<std::fabs(z) <<"  x_range_="<<x_range_<<std::endl;
	if ( std::fabs(z) > x_range_ ) {
		return well_depth_;
	}
	return 0.0;
}

Real
SquareWell2Func::dfunc( Real const /*x*/ ) const
{
	return 0.0; //This is bad news for the minimizer...
}

void
SquareWell2Func::read_data( std::istream& in ) {
	in >> x0_ >> x_range_ >>well_depth_;
	// check for optional DEGREES tag
	while ( in.peek() == ' ' || in.peek() == '\t' ) { in.get(); }
	if ( in.peek() == 'D' ) {
		std::string degrees;
		in >> degrees;
		if ( degrees == "DEGREES" ) {
			x0_      *= numeric::constants::r::deg2rad;
			x_range_ *= numeric::constants::r::deg2rad;
		}
	}
}

void
SquareWell2Func::show_definition( std::ostream &out ) const {
	out << "SQUARE_WELL_2 " << x0_ << " " << x_range_ << " " << well_depth_ << std::endl;
}

Size
SquareWell2Func::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "SQUARE_WELL_2 " <<  ( x < x0_ ) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::SquareWell2Func::SquareWell2Func() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SquareWell2Func::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( x_range_ ) ); // Real
	arc( CEREAL_NVP( well_depth_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SquareWell2Func::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( x_range_ ); // Real
	arc( well_depth_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SquareWell2Func );
CEREAL_REGISTER_TYPE( core::scoring::func::SquareWell2Func )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SquareWell2Func )
#endif // SERIALIZATION
