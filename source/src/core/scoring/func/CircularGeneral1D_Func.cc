// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
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

namespace core {
namespace scoring {
namespace constraints {

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

	while( getline( stream, line ) ){

		std::istringstream l( line );
		l >> x >> val;

		count++;
		if ( count == 1 ) xmin_ = x;
		if ( count == 2 ) xbin_ = x - xmin_;
		//if ( count > 2 ) assert( (x - x_prev) == xbin_ );

		all_vals.push_back( val );

		//x_prev = x;  // set but never used ~Labonte
	}

	num_bins_ = all_vals.size();
	data_.dimension( num_bins_ );

	for ( Size i = 1; i <= num_bins_; i++ ) data_( i ) = all_vals[ i ];
	//		std::cout << "READ: " << num_bins_ << " from " << filename << "   --> " <<
	//			" " << xmin_ << " " << xbin_ << " " << all_vals[ num_bins_ ] << std::endl;

	stream.close();
}

Real
CircularGeneral1D_Func::func( Real const x ) const {

	Real bin_real =  ( x - xmin_ ) / xbin_;

	Real const bin_wrap_real = bin_real - num_bins_ * floor( bin_real / num_bins_ ) + 1;

	assert( bin_wrap_real >= 1 && bin_wrap_real < num_bins_+1 );

	Size const bin = static_cast< Size >( bin_wrap_real );
	Real const leftover = bin_wrap_real - bin;

	Size next_bin = bin + 1;
	if ( next_bin > num_bins_ ) next_bin = 1; //wrap around.

	//	runtime_assert( bin >= 1      && bin <= data_.size() );
	//	runtime_assert( next_bin >= 1 && next_bin <= data_.size() );

	return  (data_( bin ) * ( 1 - leftover ))   +   (data_( next_bin ) * leftover) ;

}

Real
CircularGeneral1D_Func::dfunc( Real const x ) const {

	Real bin_real =  ( x - xmin_ ) / xbin_;

	Real const bin_wrap_real = bin_real - num_bins_ * floor( bin_real / num_bins_ ) + 1;
	assert( bin_wrap_real >= 1 && bin_wrap_real < num_bins_ + 1 );

	Size const bin = static_cast< Size >( bin_wrap_real );
	//Real const leftover = bin_wrap_real - bin;

	Size next_bin = bin + 1;
	if ( next_bin > num_bins_ ) next_bin = 1; //wrap around.

	return  ( data_( next_bin ) - data_( bin ) )/ xbin_;

}

} // namespace constraints
} // namespace scoring
} // namespace core
