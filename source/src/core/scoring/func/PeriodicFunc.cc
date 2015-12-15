// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/PeriodicFunc.cc
/// @brief Definition for periodic functions
/// @author Florian Richter, floric@u.washington.edu

#include <core/scoring/func/PeriodicFunc.hh>
#include <math.h>
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


bool PeriodicFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	PeriodicFunc const & other_downcast( static_cast< PeriodicFunc const & > (other) );
	if ( x0_         != other_downcast.x0_         ) return false;
	if ( k_          != other_downcast.k_          ) return false;
	if ( n_periodic_ != other_downcast.n_periodic_ ) return false;
	if ( C_          != other_downcast.C_          ) return false;

	return true;
}

bool PeriodicFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< PeriodicFunc const * > ( &other );
}

Real
PeriodicFunc::func( Real const x ) const
{
	return ( k_ * (cos( n_periodic_ * ( x - x0_ ) ) ) ) + C_;
}

Real
PeriodicFunc::dfunc( Real const x ) const
{
	return -1. * k_ * n_periodic_ * sin( n_periodic_ * (x - x0_ ) );
}

void
PeriodicFunc::read_data( std::istream& in )
{
	in >> x0_ >> n_periodic_ >> k_ >> C_;
}

void
PeriodicFunc::show_definition(std::ostream &out ) const {
	out << "PERIODIC " << x0_ << " " << n_periodic_ << " " << k_ << " " << C_
		<< std::endl;
}

//copied from HarmonicFunc.cc
Size
PeriodicFunc::show_violations(
	std::ostream& out, Real x, Size verbose_level, Real threshold
) const {
	if ( verbose_level > 100 ) {
		out << "PERIODIC " <<  func(x) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);

}


} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::PeriodicFunc::PeriodicFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::PeriodicFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( x0_ ) ); // Real
	arc( CEREAL_NVP( k_ ) ); // Real
	arc( CEREAL_NVP( n_periodic_ ) ); // Real
	arc( CEREAL_NVP( C_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::PeriodicFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( x0_ ); // Real
	arc( k_ ); // Real
	arc( n_periodic_ ); // Real
	arc( C_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::PeriodicFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::PeriodicFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_PeriodicFunc )
#endif // SERIALIZATION
