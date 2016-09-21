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
#include <numeric/MathNTensor.fwd.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathNTensorBase.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// C++ headers
#include <math.h>
#include <iostream>
#include <memory>

namespace numeric {

template< class T, numeric::Size N >
class MathNTensor : public MathNTensorBase<T>
{
public:
	typedef numeric::Size Size;
	typedef MathNTensorBase<T> parent;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	MathNTensor() :
		parent()
	{
		size_ = 0;
		data_ = nullptr;
		//std::fill( data_, data_ + size_, 0 );
	}

	//! construct Tensor from vector of dims and single element
	MathNTensor(
		utility::fixedsizearray1< Size, N > const & n_bins,
		T    const & value = T()
	) :
		parent(N),
		n_bins_( n_bins )
	{
		size_ = 1;
		for ( Size i = 1; i <= N; ++i ) {
			size_ *= n_bins[ i ];
		}
		if ( size_ == 0 ) data_ = 0;
		else data_ = new T[ size_ ];
		std::fill( data_, data_ + size_, value );
	}

	MathNTensor( MathNTensor const & src ) :
		parent( src.num_dimensions() ),
		n_bins_( src.n_bins_ ),
		size_( src.size_ ),
		data_( new T[ size_ ] )
	{
		std::copy( src.data_, src.data_ + size_, data_ ); //for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = src.data_[ ii ]; }
	}

	MathNTensor< T, N > &
	operator= ( MathNTensor< T, N > const & rhs ) {
		if ( this != &rhs ) {
			if ( size_ != rhs.size_ ) {
				delete[] data_;
				parent::set_dimensionality( rhs.num_dimensions() );
				n_bins_ = rhs.n_bins_;
				size_ = rhs.size_;
				data_ = new T[size_];
			}
			std::copy( rhs.data_, rhs.data_ + rhs.size_, data_ );
		}
		return *this;
	}


	~MathNTensor() {
		delete[] data_;
	}

	T sum() const
	{
		return std::accumulate( data_, data_ + size_, 0 );
	}

	inline
	bool is_valid_position( utility::fixedsizearray1< Size, N > const & positions ) const {
		bool validity = true;
		for ( Size i = 1; i <= N; ++i ) {
			validity = validity && ( positions[ i ] < n_bins_[ i ] );
		}
		return validity;
	}

	inline
	bool is_valid_position( utility::vector1< Size > const & positions ) const {
		bool validity = true;
		for ( Size i = 1; i <= N; ++i ) {
			validity = validity && ( positions[ i ] < n_bins_[ i ] );
		}
		return validity;
	}

	// Raw pointer constructor.  Avoid this.
	//template< Size N >
	MathNTensor(
		utility::fixedsizearray1< Size, N > const & n_bins,
		T const * data
	) :
		parent(N),
		n_bins_( n_bins )
	{
		size_ = 1;
		for ( Size i = 1; i <= N; ++i ) {
			size_ *= n_bins_[ i ];
		}
		data_ = new T[ size_ ];
		for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = data[ ii ]; }
	}

	/// @brief return number of bins for ith dimension
	inline
	Size n_bins( Size i ) const
	{
		return n_bins_[ i ];
	}

	/// @brief return number of dimensions overall
	//template< Size N >
	inline
	Size num_dimensions() const
	{
		return N;
	}

	/// @brief copies elements of argument matrix into this object at position ( layer )
	/// @details layer is the FIRST index.
	//template< Size N >
	void
	replace_layer( Size layer, MathNTensor< T, N-1 > const & matrix) {

		utility::fixedsizearray1< Size, N > position( 0 );
		position[ 1 ] = layer;
		assert( is_valid_position( position ) );

		for ( Size i = 1; i <= N - 1; ++i ) {
			assert( matrix.n_bins( i ) == n_bins_[ i + 1 ] );
		}

		// The layer offset is how many elements need to be skipped to get to
		// the beginning of the target layer.
		Size layeroffset = layer;
		for ( Size i = 1; i <= N - 1; ++i ) {
			layeroffset *= matrix.n_bins( i );
		}

		Size p = 2;
		utility::vector1< Size > indices( N+1, 0 );
		utility::fixedsizearray1< Size, N+1 > loopmaxes;
		for ( Size ii = 1; ii <= N; ++ii ) loopmaxes[ ii ] = n_bins_[ ii ];

		while ( indices[ N+1 ] == 0 ) {
			Size offset = layeroffset;
			Size slice = 1;

			for ( Size i = N; i >= 2; --i ) {
				offset += ( indices[ i ] ) * slice; // subtracting 1 bc position is 1-max not 0-max-1
				slice *= n_bins_[ i ];
			}

			utility::vector1< Size > matindices;
			for ( Size i = 2; i <= N; ++i ) {
				matindices.push_back( indices[i] );
			}
			//std::cout << offset << " " << matindices[1] << " " << matindices[2] << " " << matindices[3] << " " << matrix( matindices ) << "\n";
			data_[ offset ] = matrix( matindices );
			//std::cout << data_[ offset ] << "\n" << std::endl;

			indices[ p ]++;
			while ( indices[ p ] == loopmaxes[ p ] ) {
				indices[ p ] = 0;
				indices[ ++p ]++;
				if ( indices[ p ] != loopmaxes[ p ] ) p = 2;
			}
		}
	}

	T const data( Size const index ) {
		assert( index < size_ );
		return data_[ index ];
	}

	void
	replace_layer( Size const layer, MathMatrix< T > const & matrix )
	{
		assert( N == 3 );
		utility::vector1<Size> position( N, 1 );
		position[ 1 ] = layer;// + 1;
		assert( is_valid_position( position ) );

		assert( matrix.get_number_rows() == n_bins_[2] && matrix.get_number_cols() == n_bins_[3] );

		// if we are calling this we know we have 3 dimensions!
		Size layeroffset = layer * n_bins_[2] * n_bins_[3];
		// copy elements
		for ( Size ii = 0; ii < n_bins_[2]; ++ii ) {
			Size layerandrowoffset = layeroffset + ii*n_bins_[3];
			for ( Size jj = 0; jj < n_bins_[3]; ++jj ) {
				data_[ layerandrowoffset + jj ] = matrix( ii, jj );
			}
		}
	}


	// Accessors, both as const and nonconst reference and via
	// fixedsizearray1, vector1, or sequences of values.
	T &
	operator() ( utility::fixedsizearray1< Size, N > const & position ) {
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = N; i >= 1; --i ) {
			index += ( position[ i ] ) * slice;
			slice *= n_bins_[ i ];
		}
		return data_[ index ];
		//col + ncols_*( row + nrows_* layer) ];
	}

	T const &
	operator() ( utility::fixedsizearray1< Size, N > const & position ) const {
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = N; i >= 1; --i ) {
			index += ( position[ i ] ) * slice;
			slice *= n_bins_[ i ];
		}
		return data_[ index ];
	}

	T &
	operator() ( utility::vector1< Size > const & position ) {
		assert( N == position.size() );
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = N; i >= 1; --i ) {
			index += ( position[ i ] ) * slice;
			slice *= n_bins_[ i ];
		}
		return data_[ index ];
	}

	T const &
	operator() ( utility::vector1< Size > const & position ) const {
		assert( N == position.size() );
		assert( is_valid_position( position ) );
		Size index = 0;
		Size slice = 1;
		// I want to maintain the same data ordering. Layer seems to be the first index and gets the most stuff multiplied to it
		for ( Size i = N; i >= 1; --i ) {
			index += ( position[ i ] ) * slice;
			slice *= n_bins_[ i ];
		}
		return data_[ index ];
	}

	T & operator()( Size const b1, Size const b2, Size const b3 )
	{
		assert( N == 3 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 ];
	}

	T const & operator()( Size const b1, Size const b2, Size const b3 ) const
	{
		assert( N == 3 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 ];
	}

	T & operator()( Size const b1, Size const b2, Size const b3, Size const b4 )
	{
		assert( N == 4 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		assert( b4 < n_bins(4) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 << " which is " << data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ] << std::endl;

		return data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ];
	}

	T const & operator()( Size const b1, Size const b2, Size const b3, Size const b4 ) const
	{
		assert( N == 4 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		assert( b4 < n_bins(4) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4  << " which is " << data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ] << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ];
	}

	void
	operator/= ( Real const divisor ) {
		for ( Size ii = 0; ii < size_; ++ii ) {
			data_[ii] /= divisor;
		}
	}

	Size size() const { return size_; }

private:
	utility::fixedsizearray1< Size, N > n_bins_;   // number of each dimension
	Size size_;    // nlayers_ * nrows_ * ncols_
	T * data_;
};


}//end namespace numeric


#endif
