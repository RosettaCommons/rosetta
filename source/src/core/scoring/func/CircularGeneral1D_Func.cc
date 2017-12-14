// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file src/core/scoring/func/CircularGeneral1D_Func.cc
/// @brief A general 1D function that can be initialized by FArray or from text file. There's stuff like this all over Rosetta.
/// @author Rhiju Das

#include <core/scoring/func/CircularGeneral1D_Func.hh>
#include <core/types.hh>

#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <iostream>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

//Constructor
CircularGeneral1D_Func::CircularGeneral1D_Func(
	ObjexxFCL::FArray1D< Real > const & data,
	Real const & xmin,
	Real const & xbin):
	data_( data ),
	xmin_( xmin ),
	xbin_( xbin ),
	num_bins_( data.size() )
{
}

//Constructor
CircularGeneral1D_Func::CircularGeneral1D_Func( std::string const & filename )
{
	//Read in data file, and fill in private data.
	utility::io::izstream stream;
	stream.open( filename );

	std::string line;
	Real x( 0.0 ), /*x_prev( 0.0 ),*/ val( 0.0 );
	Size count( 0 );
	utility::vector1< Real > all_vals;

	while ( getline( stream, line ) ) {

		std::istringstream l( line );
		l >> x >> val;

		count++;
		if ( count == 1 ) xmin_ = x;
		if ( count == 2 ) xbin_ = x - xmin_;
		all_vals.push_back( val );
	}

	num_bins_ = all_vals.size();
	data_.dimension( num_bins_ );

	for ( Size i = 1; i <= num_bins_; i++ ) data_( i ) = all_vals[ i ];

	stream.close();
}

bool CircularGeneral1D_Func::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< CircularGeneral1D_Func const & > (other) );
	if ( data_     != other_downcast.data_     ) return false;
	if ( xmin_     != other_downcast.xmin_     ) return false;
	if ( xbin_     != other_downcast.xbin_     ) return false;
	if ( num_bins_ != other_downcast.num_bins_ ) return false;

	return true;
}

bool CircularGeneral1D_Func::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< CircularGeneral1D_Func const * > ( &other );
}

Real
CircularGeneral1D_Func::func( Real const x ) const {

	Real bin_real =  ( x - xmin_ ) / xbin_;
	Real const bin_wrap_real = bin_real - num_bins_ * floor( bin_real / num_bins_ ) + 1;
	debug_assert( bin_wrap_real >= 1 && bin_wrap_real < num_bins_+1 );

	auto const bin = static_cast< Size >( bin_wrap_real );
	Real const leftover = bin_wrap_real - bin;

	Size next_bin = bin + 1;
	if ( next_bin > num_bins_ ) next_bin = 1; //wrap around.

	return  (data_( bin ) * ( 1 - leftover ))   +   (data_( next_bin ) * leftover) ;
}

Real
CircularGeneral1D_Func::dfunc( Real const x ) const {
	Real bin_real =  ( x - xmin_ ) / xbin_;

	Real const bin_wrap_real = bin_real - num_bins_ * floor( bin_real / num_bins_ ) + 1;
	debug_assert( bin_wrap_real >= 1 && bin_wrap_real < num_bins_ + 1 );

	auto const bin = static_cast< Size >( bin_wrap_real );

	Size next_bin = bin + 1;
	if ( next_bin > num_bins_ ) next_bin = 1; //wrap around.

	return  ( data_( next_bin ) - data_( bin ) )/ xbin_;
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::CircularGeneral1D_Func::CircularGeneral1D_Func() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::CircularGeneral1D_Func::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( data_ ) ); // ObjexxFCL::FArray1D<core::Real>
	arc( CEREAL_NVP( xmin_ ) ); // Real
	arc( CEREAL_NVP( xbin_ ) ); // Real
	arc( CEREAL_NVP( num_bins_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::CircularGeneral1D_Func::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( data_ ); // ObjexxFCL::FArray1D<core::Real>
	arc( xmin_ ); // Real
	arc( xbin_ ); // Real
	arc( num_bins_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::CircularGeneral1D_Func );
CEREAL_REGISTER_TYPE( core::scoring::func::CircularGeneral1D_Func )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_CircularGeneral1D_Func )
#endif // SERIALIZATION
