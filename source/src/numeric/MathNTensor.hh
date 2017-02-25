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
/// *****NOTE**** The MathNTensor class is indexed at 0!
///
/// Check out MathNTensor_io.hh for file i/o functions.
///
/// @references
/// Nils Woetzl
/// Jens Meiler
///
/// @author Steven Combs, Nils Woetzl, Jens Meiler
/// @author ported to Rosetta by Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author generalized to N dimensions by Andrew Watkins
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added default constructor
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

	/// @brief Default constructor -- make a MathNTensor that's 0-dimensional with 0 elements.
	///
	MathNTensor() :
		parent()
	{
		size_ = 0;
		data_ = nullptr;
	}

	/// @brief Constructor: make a MathNTensor from fixedsizearray of dims and single element.
	/// @details The single element gets copied to every entry in the MathNTensor, initializing it.
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

	/// @brief Constructor: make a MathNTensor from vector of dims and single element.
	/// @details The single element gets copied to every entry in the MathNTensor, initializing it.
	MathNTensor(
		utility::vector1< Size > const & n_bins,
		T    const & value = T()
	) :
		parent(N),
		n_bins_( )
	{
		runtime_assert( n_bins.size() == N );
		utility::fixedsizearray1< Size, N > nbinsarray;

		size_ = 1;
		for ( Size i = 1; i <= N; ++i ) {
			nbinsarray[i] = n_bins[i];
			size_ *= n_bins[ i ];
		}
		n_bins_ = nbinsarray;
		if ( size_ == 0 ) data_ = 0;
		else data_ = new T[ size_ ];
		std::fill( data_, data_ + size_, value );
	}

	/// @brief Copy constructor.
	///
	MathNTensor( MathNTensor const & src ) :
		parent( src.num_dimensions() ),
		n_bins_( src.n_bins_ ),
		size_( src.size_ ),
		data_( new T[ size_ ] )
	{
		std::copy( src.data_, src.data_ + size_, data_ ); //for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = src.data_[ ii ]; }
	}

	/// @brief Assignment operator.
	///
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

	/// @brief Destructor.
	///
	~MathNTensor() {
		delete[] data_;
	}

	T sum() const
	{
		return std::accumulate( data_, data_ + size_, 0 );
	}

	/// @brief Is the position in the N-Tensor given by the positions vector a valid position?
	/// @details Returns true if each coordinate is in the range [1, dimension], false otherwise.
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

	/// @brief Raw pointer constructor.  Avoid this.
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

	/// @brief Return number of dimensions overall (i.e. the dimensionality
	/// of this MathNTensor).
	inline
	Size num_dimensions() const
	{
		return N;
	}

	/// @brief Copies elements of argument matrix into this object at position ( layer )
	/// @details layer is the FIRST index.
	void
	replace_layer( Size layer, MathNTensor< T, N-1 > const & matrix) {

		utility::fixedsizearray1< Size, N > position( 0 );
		position[ 1 ] = layer;
		runtime_assert( is_valid_position( position ) );

		for ( Size i = 1; i <= N - 1; ++i ) {
			runtime_assert( matrix.n_bins( i ) == n_bins_[ i + 1 ] );
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

	T const & data( Size const index ) const {
		runtime_assert( index < size_ );
		return data_[ index ];
	}

	/// @brief ONLY for 3-dimensional tensors, replace a layer with a 2D matrix.
	///
	void
	replace_layer( Size const layer, MathMatrix< T > const & matrix )
	{
		runtime_assert( N == 3 );
		utility::vector1<Size> position( N, 1 );
		position[ 1 ] = layer;// + 1;
		runtime_assert( is_valid_position( position ) );

		runtime_assert( matrix.get_number_rows() == n_bins_[2] && matrix.get_number_cols() == n_bins_[3] );

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

	/// @brief ONLY for 2-dimensional tensors, give me a MathMatrix that's identical to the tensor.  (That is, convert
	/// the MathNTensor to MathMatrix form.
	/// @details Throws if this is not a 2D tensor.
	/// @return A MathMatrix containing the full contents of the MathNTensor.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	MathMatrix< T >
	get_mathmatrix()
	{
		runtime_assert_string_msg( N == 2, "Error in numeric::MathNTensor::get_MathMatrix(): This operation can only be used for two-dimensional tensors." );
		MathMatrix < T > outmatrix( n_bins_[1], n_bins_[2] );
		Size count(0);
		for ( Size i=0, imax=n_bins_[1]; i<imax; ++i ) {
			for ( Size j=0, jmax=n_bins_[2]; j<jmax; ++j ) {
				outmatrix(i,j)=data_[count];
				++count;
			}
		}
		return outmatrix;
	}


	/// @brief Access a position in the N-dimensional tensor.
	/// @details Accessors, both as const and nonconst reference and via
	/// fixedsizearray1, vector1, or sequences of values.
	T &
	operator() ( utility::fixedsizearray1< Size, N > const & position ) {
		runtime_assert_string_msg( is_valid_position( position ), "Error in numeric::MathNTensor::operator(): Attempting to access an invalid entry in the tensor." );
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

	/// @brief Use an N-vector of coordinates to access a position in the N-dimensional tensor.
	/// @details Const access (assignment prohibited).
	T const &
	operator() ( utility::fixedsizearray1< Size, N > const & position ) const {
		runtime_assert_string_msg( is_valid_position( position ), "Error in numeric::MathNTensor::operator(): Attempting to access an invalid entry in the tensor." );
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
		runtime_assert_string_msg( N == position.size(), "Error in numeric::MathNTensor::operator(): Attempting to access a tensor entry using an array of invalid size." );
		runtime_assert_string_msg( is_valid_position( position ), "Error in numeric::MathNTensor::operator(): Attempting to access an invalid entry in the tensor." );
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
		runtime_assert_string_msg( N == position.size(), "Error in numeric::MathNTensor::operator(): Attempting to access a tensor entry using an array of invalid size." );
		runtime_assert_string_msg( is_valid_position( position ), "Error in numeric::MathNTensor::operator(): Attempting to access an invalid entry in the tensor." );
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
		runtime_assert( N == 3 );
		runtime_assert( b1 < n_bins(1) );
		runtime_assert( b2 < n_bins(2) );
		runtime_assert( b3 < n_bins(3) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 ];
	}

	T const & operator()( Size const b1, Size const b2, Size const b3 ) const
	{
		runtime_assert( N == 3 );
		runtime_assert( b1 < n_bins(1) );
		runtime_assert( b2 < n_bins(2) );
		runtime_assert( b3 < n_bins(3) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) + b2 * n_bins(3) + b3 ];
	}

	T & operator()( Size const b1, Size const b2, Size const b3, Size const b4 )
	{
		runtime_assert( N == 4 );
		runtime_assert( b1 < n_bins(1) );
		runtime_assert( b2 < n_bins(2) );
		runtime_assert( b3 < n_bins(3) );
		runtime_assert( b4 < n_bins(4) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 << " which is " << data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ] << std::endl;

		return data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ];
	}

	T const & operator()( Size const b1, Size const b2, Size const b3, Size const b4 ) const
	{
		runtime_assert( N == 4 );
		runtime_assert( b1 < n_bins(1) );
		runtime_assert( b2 < n_bins(2) );
		runtime_assert( b3 < n_bins(3) );
		runtime_assert( b4 < n_bins(4) );
		//std::cout << "Fetching element " << b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4  << " which is " << data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ] << std::endl;
		return data_[ b1 * n_bins(2) * n_bins(3) * n_bins(4) + b2 * n_bins(3) * n_bins(4) + b3 * n_bins(4) + b4 ];
	}

	T const & operator()( Size const b1, Size const b2, Size const b3, Size const b4, Size const b5 ) const
	{
		assert( N == 5 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		assert( b4 < n_bins(4) );
		assert( b5 < n_bins(5) );
		utility::fixedsizearray1< Size, N > position;
		position[ 1 ] = b1;
		position[ 2 ] = b2;
		position[ 3 ] = b3;
		position[ 4 ] = b4;
		position[ 5 ] = b5;
		return (*this)( position );
	}


	T const & operator()( Size const b1, Size const b2, Size const b3, Size const b4, Size const b5, Size const b6 ) const
	{
		assert( N == 6 );
		assert( b1 < n_bins(1) );
		assert( b2 < n_bins(2) );
		assert( b3 < n_bins(3) );
		assert( b4 < n_bins(4) );
		assert( b5 < n_bins(5) );
		assert( b6 < n_bins(6) );
		utility::fixedsizearray1< Size, N > position;
		position[ 1 ] = b1;
		position[ 2 ] = b2;
		position[ 3 ] = b3;
		position[ 4 ] = b4;
		position[ 5 ] = b5;
		position[ 6 ] = b6;
		return (*this)( position );
	}

	void
	operator/= ( Real const divisor ) {
		for ( Size ii = 0; ii < size_; ++ii ) {
			data_[ii] /= divisor;
		}
	}

	Size size() const { return size_; }

	utility::fixedsizearray1< Size, N > const & n_bins() const { return n_bins_; }

private:
	utility::fixedsizearray1< Size, N > n_bins_;   // number of each dimension

	/// @brief The total number of bins in the MathNTensor (the product of all dimensions).
	///
	Size size_;    // nlayers_ * nrows_ * ncols_

	/// @brief Array of data stored in this N-tensor.
	///
	T * data_;
};

}//end namespace numeric


#endif
