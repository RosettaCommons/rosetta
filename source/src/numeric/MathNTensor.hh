// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// construction/destructor of 3-D Matrix's with some functions
///
/// @details
/// This is an implementation of an algorithm that was taken from BCL (Jens Meiler)
/// *****NOTE**** The MathTensor class is indexed at 0!!!!
///
/// @references
/// Nils Woetzl
/// Jens Meiler
///
/// @author Steven Combs, Nils Woetzl, Jens Meiler
/// @author ported to Rosetta by Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author generalized to N dimensions by Andrew Watkins
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathNTensor_hh
#define INCLUDED_numeric_MathNTensor_hh

// Package headers
#include <numeric/types.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ headers
#include <math.h>
#include <iostream>

namespace numeric {

template< typename T >
class MathNTensor
{
public:


	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////


	//! construct Tensor from vector of dims and single element
	MathNTensor(
		utility::vector1< Size > const & n_dimensions,
		T    const & value = T()
	) :
		n_dimensions_( n_dimensions )
	{
		n_xs_ = n_dimensions_.size();
		size_ = 1;
		for ( Size i = 1; i <= n_dimensions.size(); ++i ) {
			size_ *= n_dimensions[ i ];
		}
		if ( size_ == 0 ) data_ = 0;
		else data_ = new T[ size_ ];
		//data_( size_ != 0 ? new T[ size_ ] : 0 );
		std::fill( data_, data_ + size_, value );
	}
	//! construct Tensor from optional number of layers, optional number of rows, optional number of columns, and optional single element
	MathNTensor( MathNTensor const & src ) :
		n_dimensions_( src.n_dimensions_ ),
		size_( src.size_ ),
		data_( new T[ size_ ] )
	{
		n_xs_ = n_dimensions_.size();
		std::copy( src.data_, src.data_ + size_, data_ ); //for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = src.data_[ ii ]; }
	}

	MathNTensor &
	operator = ( MathNTensor const & rhs ) {
		if ( this != & rhs ) {
			if ( size_ != rhs.size_ ) {
				delete[] data_;
			}
			n_dimensions_ = rhs.n_dimensions_;
			n_xs_ = n_dimensions_.size();
			size_ = rhs.size_;
			data_ = new T[size_];
			std::copy( rhs.data_, rhs.data_ + rhs.size_, data_ );
		}
		return *this;
	}


	~MathNTensor() {
		delete [] data_;
	}

	inline
	bool is_valid_position( utility::vector1< Size > positions ) const {
		bool validity = true;
		for ( Size i = 1; i <= n_dimensions_.size(); ++i ) {
			validity = validity && ( positions[ i ] <= n_dimensions_[ i ] && positions[ i ] >= 1 );
		}
		return validity;
	}


	// Raw pointer constructor.  Avoid this.
	MathNTensor(
		utility::vector1< Size > const & dims,
		T const * data
	) :
		n_dimensions_( dims )
	{
		size_ = 1;
		for ( Size i = 1; i <= n_dimensions_.size(); ++i ) {
			size_ *= n_dimensions_[ i ];
		}
		data_ = new T[ size_ ];
		for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = data[ ii ]; }
		n_xs_ = n_dimensions_.size();
	}

	/// @brief return number of bins for ith dimension
	Size n_dimensions( Size i ) const
	{
		return n_dimensions_[ i ];
	}

	/// @brief return number of dimensions overall
	Size num_dimensions() const
	{
		return n_xs_;
	}

	/// @brief copies elements of argument matrix into this object at position ( layer )
	void
	replace_layer( Size layer, MathNTensor< T > const & matrix) {

		utility::vector1<Size> position( n_xs_, 1 );
		position[ 1 ] = layer + 1;
		assert( is_valid_position( position ) );

		for ( Size i = 1; i <= n_xs_-1; ++i ) {
			assert( matrix.n_dimensions( i ) == n_dimensions_[ i + 1 ] );
		}

		Size layeroffset = layer;
		for ( Size i = 1; i <= n_xs_-1; ++i ) {
			layeroffset *= matrix.n_dimensions( i );
		}

		Size p = 2;
		utility::vector1< Size > indices( n_xs_+1, 0 );
		while ( indices[ n_xs_ ] == 0 ) {
			Size thisoffset = layeroffset;

			utility::vector1< Size > matindices;
			for ( Size i = 2; i <= n_xs_; ++i ) {
				thisoffset += indices[i] * n_dimensions_[i];
				matindices.push_back( indices[i]+1 );
			}
			data_[ thisoffset + indices[ n_xs_ ] ] = matrix( matindices );

			indices[ p ]++;
			while ( indices[ n_xs_ ] == n_dimensions_[ n_xs_ ] ) {
				indices[ p ] = 0;
				indices[ ++p ]++;
				if ( indices[ n_xs_ ] != n_dimensions_[ n_xs_ ] ) p = 2;
			}
		}

		// copy elements
		//for ( Size ii = 0; ii < nrows_; ++ii ) {
		// Size layerandrowoffset = layeroffset + ii*ncols_;
		// for ( Size jj = 0; jj < ncols_; ++jj ) {
		//  data_[ layerandrowoffset + jj ] = matrix( ii, jj );
		// }
		//}
	}

	void
	replace_layer( Size const layer, MathMatrix< T > const & matrix )
	{
		utility::vector1<Size> position( n_xs_, 1 );
		position[ 1 ] = layer + 1;
		assert( is_valid_position( position ) );

		assert( matrix.get_number_rows() == n_dimensions_[2] && matrix.get_number_cols() == n_dimensions_[3] );

		// if we are calling this we know we have 3 dimensions!
		Size layeroffset = layer * n_dimensions_[2] * n_dimensions_[3];
		// copy elements
		for ( Size ii = 0; ii < n_dimensions_[2]; ++ii ) {
			Size layerandrowoffset = layeroffset + ii*n_dimensions_[3];
			for ( Size jj = 0; jj < n_dimensions_[3]; ++jj ) {
				data_[ layerandrowoffset + jj ] = matrix( ii, jj );
			}
		}
	}


	T &
	operator() ( utility::vector1< Size > position ) {
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = n_xs_; i >= 1; --i ) {
			index += ( position[ i ]-1 ) * slice; // subtracting 1 bc position is 1-max not 0-max-1
			slice *= n_dimensions_[ i ];
		}
		return data_[ index ];
		//col + ncols_*( row + nrows_* layer) ];
	}

	T const &
	operator() ( utility::vector1< Size > position ) const {
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = n_xs_; i >= 1; --i ) {
			index += ( position[ i ]-1 ) * slice; // subtracting 1 bc position is 1-max not 0-max-1
			slice *= n_dimensions_[ i ];
		}
		return data_[ index ];
	}

private:
	Size n_xs_;
	utility::vector1< Size > n_dimensions_;   // number of each dimension
	Size size_;    // nlayers_ * nrows_ * ncols_
	T * data_;

};


}//end namespace numeric


#endif
