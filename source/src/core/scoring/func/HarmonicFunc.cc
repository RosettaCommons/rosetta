// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/HarmonicFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/func/HarmonicFunc.hh>

#include <core/types.hh>



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

FuncOP
HarmonicFunc::clone() const
{
	return utility::pointer::make_shared< HarmonicFunc >( x0_, sd_ );
}

bool HarmonicFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< HarmonicFunc const & > (other) );
	if ( x0_ != other_downcast.x0_ ) return false;
	if ( sd_ != other_downcast.sd_ ) return false;

	return true;
}

bool HarmonicFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< HarmonicFunc const * > ( &other ) != nullptr;
}


Real
HarmonicFunc::func( Real const x ) const
{
	Real const z = ( x-x0_ )/sd_;
	return z * z;
}

Real
HarmonicFunc::dfunc( Real const x ) const
{
	return 2 * (x-x0_) / ( sd_ * sd_ );
}

void
HarmonicFunc::read_data( std::istream& in ) {
	in >> x0_ >> sd_;
}

void
HarmonicFunc::show_definition( std::ostream &out ) const {
	out << "HARMONIC " << x0_ << " " << sd_ << std::endl;
}

Size
HarmonicFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if ( verbose_level > 100 ) {
		out << "HARM " <<  ( x - x0_ )/sd_ << std::endl;
	} else if ( verbose_level > 70 ) {
		if ( x < x0_  && ( this->func(x) > threshold ) ) out << "-";
		else if ( x > x0_ && ( this->func(x) > threshold ) ) out << "+";
		else out << ".";
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::HarmonicFunc::HarmonicFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::HarmonicFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( sd_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::HarmonicFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( sd_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::HarmonicFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::HarmonicFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_HarmonicFunc )
#endif // SERIALIZATION
