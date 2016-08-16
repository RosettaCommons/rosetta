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
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathTensor_hh
#define INCLUDED_numeric_MathTensor_hh

// Package headers
#include <numeric/types.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <math.h>
#include <iostream>

namespace numeric {

template< typename T >
class MathTensor
{
public:


	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////


	//! construct Tensor from optional number of layers, optional number of rows, optional number of columns, and optional single element
	MathTensor(
		Size const layers = 0,
		Size const rows = 0,
		Size const cols = 0,
		T    const & value = T()
	) :
		nlayers_( layers ),
		nrows_( rows ),
		ncols_( cols ),
		size_( nlayers_ * nrows_ * ncols_ ),
		data_( size_ != 0 ? new T[ size_ ] : 0 )
	{
		std::fill( data_, data_ + size_, value );
	}
	//! construct Tensor from optional number of layers, optional number of rows, optional number of columns, and optional single element
	MathTensor( MathTensor const & src ) :
		nlayers_( src.nlayers_ ),
		nrows_( src.nrows_ ),
		ncols_( src.ncols_ ),
		size_( src.size_ ),
		data_( new T[ size_ ] )
	{
		std::copy( src.data_, src.data_ + size_, data_ ); //for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = src.data_[ ii ]; }
	}

	MathTensor &
	operator = ( MathTensor const & rhs ) {
		if ( this != & rhs ) {
			if ( size_ != rhs.size_ ) {
				delete[] data_;
			}
			nlayers_ = rhs.nlayers_;
			ncols_ = rhs.ncols_;
			nrows_ = rhs.nrows_;
			size_ = rhs.size_;
			data_ = new T[size_];
			std::copy( rhs.data_, rhs.data_ + rhs.size_, data_ );
		}
		return *this;
	}


	~MathTensor() {
		delete [] data_;
	}

	inline
	bool is_valid_position( Size layer, Size row, Size column ) const {
		return ( layer < nlayers_ )
			&& ( row   < nrows_   )
			&& ( column< ncols_   );
	}


	// Raw pointer constructor.  Avoid this.
	MathTensor(
		Size const layers,
		Size const rows,
		Size const cols,
		T const * data
	) :
		nlayers_( layers ),
		nrows_( rows ),
		ncols_( cols ),
		size_( layers * rows * cols ),
		data_( new T[ size_ ] )
	{
		for ( Size ii = 0; ii < size_; ++ii ) { data_[ ii ] = data[ ii ]; }
	}

	/// @brief return number of layers
	Size nlayers() const
	{
		return nlayers_;
	}

	/// @brief return number of rows
	Size nrows() const
	{
		return nrows_;
	}

	//! return number of columns
	Size ncols() const
	{
		return ncols_;
	}

	/// @brief copies elements of argument matrix into this object at position ( layer )
	void
	replace_layer( Size const layer, MathMatrix< T > const & matrix) {
		assert( is_valid_position( layer, 0, 0));
		assert( matrix.get_number_rows() == nrows_ && matrix.get_number_cols() == ncols_ );

		Size layeroffset = layer * nrows_ * ncols_;
		// copy elements
		for ( Size ii = 0; ii < nrows_; ++ii ) {
			Size layerandrowoffset = layeroffset + ii*ncols_;
			for ( Size jj = 0; jj < ncols_; ++jj ) {
				data_[ layerandrowoffset + jj ] = matrix( ii, jj );
			}
		}

	}

	T &
	operator() ( Size layer, Size row, Size col ) {
		assert( is_valid_position( layer, row, col ) );
		return data_[ col + ncols_*( row + nrows_* layer) ];
	}

	T const &
	operator() ( Size layer, Size row, Size col ) const {
		assert( is_valid_position( layer, row, col ) );
		return data_[ col + ncols_*( row + nrows_* layer) ];
	}

private:

	Size nlayers_; // number of layers
	Size nrows_;   // number of rows
	Size ncols_;   // number columns
	Size size_;    // nlayers_ * nrows_ * ncols_
	T * data_;

};


}//end namespace numeric


#endif

