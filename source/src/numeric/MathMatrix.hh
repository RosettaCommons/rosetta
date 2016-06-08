// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// construction/destructor of Matrix's with some functions
///
/// @details
/// This is an implementation of an algorithm that was taken from BCL (Jens Meiler)
/// The Matrix is constructed out of arrays and places values into rows/columns based on
/// however many columns/rows you specify. Actual operations of the MathMatrix are implemented
/// in numeric/MathMatrix_operations.hh. To access specific values (elements), you must use
/// the operator (). For example: to access row 5, column 3 of a matrix, you would use
/// matrix(5,3). *****NOTE**** The MathMatrix class is indexed at 0!!!!
///
/// @references
/// Nils Woetzl
/// Jens Meiler
///
/// @author Steven Combs, Nils Woetzl, Jens Meiler
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathMatrix_hh
#define INCLUDED_numeric_MathMatrix_hh

// Package headers
#include <numeric/types.hh>
#include <numeric/MathVector.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <math.h>
#include <iostream>
#include <utility/assert.hh>

namespace numeric {
template<typename T>
class MathMatrix
{
public:

	//////////
	// data //
	//////////


	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief default constructor
	MathMatrix< T>() :
	NumberRows_( 0),
	NumberCols_( 0),
	size_( 0 ),
	data_( NULL)
	{
	}

	/// @brief construct from dimension and possible filler
	/// @param NUMBER_ROWS number of rows in matrix
	/// @param NUMBER_COLS number of cols in matrix
	/// @param FILL_VALUE assign every element to that value
	explicit MathMatrix< T>
	(
		const Size NUMBER_ROWS,
		const Size NUMBER_COLS,
		const T &FILL_VALUE = T( 0)
	) :
	NumberRows_( NUMBER_ROWS),
	NumberCols_( NUMBER_COLS),
	size_( NumberRows_ * NumberCols_ ),
	data_( new T[ size_ ])
	{
		// set all values to FILL_VALUE
		std::fill( data_, data_ + size_, FILL_VALUE);
	}

	/// @brief construct from dimension and pointer to data
	/// @param NUMBER_ROWS number of rows in matrix
	/// @param NUMBER_COLS number of cols in matrix
	/// @param DATA pointer to field of data
	MathMatrix< T>
	(
		const Size NUMBER_ROWS,
		const Size NUMBER_COLS,
		const T *DATA
	) :
		NumberRows_( NUMBER_ROWS),
		NumberCols_( NUMBER_COLS),
		size_( NumberRows_ * NumberCols_ ),
		data_( new T[ NumberRows_ * NumberCols_])
	{
		// copy data
		std::copy( DATA, DATA + size_, data_);
	}

	/// @brief copy constructor from Matrix
	/// @param MATRIX matrix to be copied from
	MathMatrix< T>( MathMatrix< T > const & MATRIX ):
		NumberRows_( MATRIX.NumberRows_ ),
		NumberCols_( MATRIX.NumberCols_ ),
		size_( MATRIX.size_ ),
		data_( new T[ size_ ] )
	{
		std::copy( MATRIX.data_, MATRIX.data_ + size_, data_);
	}


	/// @brief Clone function
	/// @return pointer to new MatrixInterface< T>
	MathMatrix< T> *Clone() const
	{
		return new MathMatrix< T>( *this);
	}

	/// @brief destructor
	~MathMatrix< T>()
	{
		delete[] data_;
	}

	/////////////////
	// data access //
	/////////////////


	/// @brief get number of rows
	/// @return number of rows
	Size get_number_rows() const
	{
		return NumberRows_;
	}

	/// @brief get number of columns
	/// @return number of columns
	Size get_number_cols() const
	{
		return NumberCols_;
	}

	/// @brief number of elements
	/// @return total number of elements in matrix
	Size get_number_elements() const
	{
		return size_;
	}

	/// @brief number of elements
	/// @return total number of elements in matrix
	Size size() const
	{
		return size_;
	}

	/// @brief pointer to First Element
	/// @return const pointer to first element in range containing all elements of Matrix
	const T *begin() const
	{
		return data_;
	}

	/// @brief pointer to First Element
	/// @return pointer to first element in range containing all elements of Matrix
	T *begin()
	{
		return data_;
	}

	/// @brief pointer to end of range
	/// @return const pointer to address one after last element in Matrix
	const T *end() const
	{
		return data_ + size_;
	}

	/// @brief pointer to end of range
	/// @return pointer to address one after last element in Matrix
	T *end()
	{
		return data_ + size_;
	}


	/// @return Row of Matrix
	MathVector< T> get_row( const Size ROW) const
	{
		return MathVector< T>( NumberCols_, operator[]( ROW));
	}

	/// @return Col of Matrix
	MathVector< T> get_col( const Size COL) const
	{
		//create a vector of the size of NumberRows
		MathVector< T> col( NumberRows_);

		//ptr to first element in Col
		T *ptr = col.begin();

		//iterate over all rows
		for ( Size i( 0); i < NumberRows_; ++i, ++ptr ) {
			( *ptr) = operator()( i, COL);
		}

		//return the column
		return col;
	}

	////////////////
	// operations //
	////////////////

	/// @brief is matrix a square matrix
	/// @return true if number of cols and rows are idnetical
	bool is_square() const
	{
		return NumberRows_ != 0 && NumberRows_ == NumberCols_;
	}

	/// @brief is matrix a diagonal matrix
	/// @return true if all but the elements in the diagonal are 0
	bool is_diagonal() const
	{
		// if matrix is not square or empty
		if ( !is_square() || get_number_elements() == 0 ) {
			return false;
		}

		// check that all but the elements in the diagonal are 0
		for ( Size i = 0; i < NumberRows_ - 1; i++ ) {
			for ( Size j = i + 1; j < NumberCols_; j++ ) {
				if ( operator()( j, i) != T( 0) || operator()( i, j) != T( 0) ) {
					return false;
				}
			}
		}

		// return true if all elements but diagonal are 0
		return true;
	}

	/// @brief is matrix a tridiagonal matrix
	/// @return if diagonal and adjecent diagonals are filled and the rest is 0
	bool is_tri_diagonal() const
	{
		// if matrix is not square and does not have at least 2 elements in each dimension
		if ( !is_square() || NumberRows_ < 2 ) {
			return false;
		}

		// check that all but the inner three diagonal elements are 0
		for ( Size i( 0); i < NumberRows_ - 2; i++ ) {
			for ( Size j( i + 2); j < NumberCols_; j++ ) {
				if ( operator()( j, i) != T( 0) || operator()( i, j) != T( 0) ) {
					return false;
				}
			}
		}

		// return true if all elements but the inner three diagonals are 0
		return true;
	}


	//////////////////////
	// Matrix functions //
	//////////////////////

	/// @brief check dimension agreement of two Matrices
	/// @param MATRIX_LHS rhs matrix
	/// @param MATRIX_RHS lhs matrix
	/// @return true if number rows and cols are the same between both Matrices
	inline bool same_dimensions
	(
		const MathMatrix< T> & MATRIX_LHS,
		const MathMatrix< T> & MATRIX_RHS
	)
	{
		return    MATRIX_LHS.get_number_rows() == MATRIX_RHS.get_number_rows()
			&& MATRIX_LHS.get_number_cols() == MATRIX_RHS.get_number_cols();
	}

	/// @brief check inverse dimension agreement of two Matrices
	/// comapre number ros of lhs with number cols of rhs and number cols of lhs with number rows of rhs
	/// @param MATRIX_LHS rhs matrix
	/// @param MATRIX_RHS lhs matrix
	/// @return true if number rows with cols and cols with rows agree between both Matrices
	inline
	bool
	inverse_dimensions
	(
		const MathMatrix< T> &MATRIX_LHS,
		const MathMatrix< T> &MATRIX_RHS
	)
	{
		return    MATRIX_LHS.get_number_rows() == MATRIX_RHS.get_number_cols()
			&& MATRIX_LHS.get_number_cols() == MATRIX_RHS.get_number_rows();
	}

	/// @brief check dimensions for multiplication A*B
	/// compare number cols of lhs with number rows of rhs
	/// @param MATRIX_LHS rhs matrix
	/// @param MATRIX_RHS lhs matrix
	/// @return true if number cols rhs and number cols lhs agree
	inline
	bool
	multiplication_dimension
	(
		const MathMatrix< T> &MATRIX_LHS,
		const MathMatrix< T> &MATRIX_RHS
	)
	{
		return MATRIX_LHS.get_number_cols() == MATRIX_RHS.get_number_rows();
	}


	MathMatrix< T> & set_zero()
	{
		// fill all with 0
		std::fill( data_, data_ + size_, T( 0));

		//end
		return *this;
	}


	/// @return Transposed of Matrix
	inline MathMatrix< T> transpose( const MathMatrix< T> &MATRIX)
	{
		return MathMatrix< T>( MATRIX).transpose();
	}


	MathMatrix< T> & inverse()
	{
		is_square() == true ? inverse_square_matrix() : inverse_rectangular_matrix();
		return *this;
	}


	/// @return transpose of matrix
	inline MathMatrix< T> & transpose()
	{
		MathMatrix< T> newthis( NumberCols_, NumberRows_);
		for ( Size i( 0); i < newthis.NumberRows_; ++i ) {
			for ( Size j( 0); j < newthis.NumberCols_; ++j ) {
				newthis( i, j) = operator()( j, i);
			}
		}

		//end
		return ( operator =( newthis));
	}

	/// @return invert rectangular matrices exactly
	inline MathMatrix< T> & inverse_rectangular_matrix()
	{
		MathMatrix< T> newmatrix( *this);
		bool transposed( NumberRows_ < NumberCols_);
		if ( !transposed ) {
			newmatrix = transpose( newmatrix) * ( newmatrix);
		} else {
			newmatrix = ( newmatrix) * ( transpose( ( newmatrix)));
		}

		newmatrix.inverse_square_matrix();

		if ( transposed ) {
			newmatrix *= *this;
			newmatrix.transpose();
		} else {
			newmatrix *= transpose( ( *this));
		}

		return ( operator =( newmatrix));
		//  return std::copy(newmatrix.data_, newmatrix.data_ + newmatrix.get_number_rows() * newmatrix.get_number_cols(), newmatrix.data_);
	}


	/// invert small square matrices exactly
	inline MathMatrix< T> inverse_square_matrix(){
		// do quick inverse if diagonal
		if ( is_diagonal() ) {
			return inverse_diagonal_matrix();
		}

		// do quick inverse if tridiagonal
		if ( is_tri_diagonal() ) {
			return inverse_tridiagonal_matrix();
		}

		MathMatrix< T> newmatrix( NumberRows_, NumberCols_);
		newmatrix.set_unit();

		for ( Size k( 0); k < NumberRows_; ++k ) {
			Size index( pivot( k));

			if ( index != 0 ) {
				newmatrix.swap_rows( k, index);
			}

			T a1 = operator()( k, k);
			for ( Size j( 0); j < NumberRows_; ++j ) {
				operator()( k, j) /= a1;
				newmatrix(  k, j) /= a1;
			}
			for ( Size i( 0); i < NumberRows_; ++i ) {
				if ( i == k ) {
					continue;
				}

				const T a2 = operator()( i, k);
				for ( Size j( 0); j < NumberRows_; ++j ) {
					operator()( i, j) -= a2 * operator()( k, j);
					newmatrix(  i, j) -= a2 * newmatrix(  k, j);
				}
			}
		}

		return ( operator =( newmatrix));
	}


	inline MathMatrix< T> & inverse_diagonal_matrix()
	{
		for ( Size i( 0); i < NumberRows_; ++i ) {
			if ( operator()( i, i) != T( 0) ) {
				operator()( i, i) = T( 1) / operator()( i, i);
			}
		}

		return *this;
	}


	/// @return this algorithm was found on this page: http://www.csit.fsu.edu/~burkardt/math2071/math2071.html
	/// invert tridiagonal matrix for all diagonal elements
	inline MathMatrix< T> & inverse_tridiagonal_matrix()
	{
		Size n( NumberRows_);
		MathMatrix< T> newmatrix( n, n);

		// LU decomposition (in lower and upper triangular matrix)
		for ( Size i( 1); i < n; ++i ) {
			operator()( i, i) -= operator()( i - 1, i) * operator()( i, i - 1) / operator()( i - 1, i - 1);
			operator()( i, i - 1) /= operator()( i - 1, i - 1);
		}

		// compute inverse form L and U
		for ( Size j( 0); j < n; ++j ) {
			// Solve L * y = b.
			MathVector< T> y( n, T( 0));
			for ( Size i( j); i < n; ++i ) {
				if ( i == j ) y( i) = T( 1);
				if ( i  > j ) y( i) -= operator()( i, i - 1) * y( i - 1);
			}

			// Solve U * x = y.
			for ( Size i( n - 1); i > 0; i-- ) {
				newmatrix( i, j) = y( i) / operator()( i, i);
				y( i - 1) -= operator()( i - 1, i) * newmatrix( i, j);
			}
			newmatrix( 0, j) = y( 0) / operator()( 0, 0);
		}

		return( operator =( newmatrix));
	}


	/// @return set all elements in matrix to T( 0) but diagonal elements to T( 1)
	MathMatrix< T> & set_unit()
	{
		//set all elelemnts to T( 0)
		set_zero();

		//iterate over diagonal of matrix and set elements to T( 1)
		for ( Size i = 0; i < NumberRows_ && i < NumberCols_; i++ ) {
			( *this)( i, i) = T( 1);
		}

		return *this;
	}

	/// @return private helper function for computing the determinante / inverting a square matrix
	inline Size pivot( const Size ROW)
	{
		Size k( ROW);

		if ( k != ROW ) {
			swap_rows( k, ROW);
			return k;
		}

		return 0;
	}

	/// @return copies elements of argument VECTOR into this object at position (ROW)
	inline MathMatrix< T> & replace_row
	(
		const Size ROW,
		const MathVector< T> &VECTOR
	)
	{
		IsValidPosition( ROW, 0);
		IsValidPosition( ROW, VECTOR.size() - 1);
		const T *dat( VECTOR.begin()), *dat_end( VECTOR.end());
		// copy elements
		for ( T *ptr( operator[]( ROW)); dat != dat_end; ++ptr, ++dat ) {

			( *ptr) = ( *dat);
		}

		return *this;
	}

	/// @return copies elements of argument VECTOR into this object at position (COL)
	inline MathMatrix< T> & replace_col
	(
		const Size COL,
		const MathVector< T> &VECTOR
	)
	{
		//check valid positions
		IsValidPosition( 0                    , COL);
		IsValidPosition( VECTOR.size() - 1 , COL);

		const T *dat( VECTOR.begin()), *dat_end( VECTOR.end());
		// copy elements
		for ( Size i( 0); i < NumberRows_ && dat != dat_end; ++i, ++dat ) {
			operator()( i, COL) = ( *dat);
		}

		//end
		return *this;
	}

	/// @return swap rows ROW_A and ROW_B
	inline MathMatrix< T> &swap_rows
	(
		const Size ROW_A,
		const Size ROW_B
	)
	{
		IsValidPosition( ROW_A, 0);
		IsValidPosition( ROW_B, 0);

		//swap each pair in rows
		for ( Size i( 0); i < NumberCols_; ++i ) {
			std::swap( operator()( ROW_A, i), operator()( ROW_B, i));
		}

		//end
		return *this;
	}

	/// @return swap columns COL_A and COL_B
	inline MathMatrix< T> & swap_cols( const Size COL_A, const Size COL_B)
	{
		IsValidPosition( 0, COL_A);
		IsValidPosition( 0, COL_B);

		//swap each pair in cols
		for ( Size i( 0); i < NumberRows_; ++i ) {
			std::swap( operator()( i, COL_A), operator()(  i, COL_B));
		}

		//end
		return *this;
	}


	//////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////


	/// @return check whether position is valid
	bool IsValidPosition( const Size ROW, const Size COL) const
	{
		if ( ROW > NumberRows_ ) {
			utility_exit_with_message("ROW extends size of Matrix");
			return false;
		}
		if ( COL > NumberCols_ ) {
			utility_exit_with_message("COL extends size of Matrix");
			return false;
		} else return true;

	}


	///////////////
	// operators //
	///////////////

	/// @brief return reference to changeable element ( ROW, COL)
	/// @param ROW the row number, starting with 0
	/// @param COL the col number, starting with 0
	/// @return changable reference to the element defined bey ROW and COL number
	T &operator()( const Size ROW, const Size COL)
	{
		assert( ROW < NumberRows_ );
		assert( COL < NumberCols_ );
		return data_[ ROW * NumberCols_ + COL];
	}

	/// @brief return reference to const  element ( ROW, COL)
	/// @param ROW the row number, starting with 0
	/// @param COL the col number, starting with 0
	/// @return const element defined bey ROW and COL number
	const T &operator()( const Size ROW, const Size COL) const
	{
		assert( ROW < NumberRows_ );
		assert( COL < NumberCols_ );
		return data_[ ROW * NumberCols_ + COL];
	}

	/// @brief assignment from Matrix
	/// @param MATRIX the matrix used as source
	/// @return reference to this Matrix
	MathMatrix< T> &operator = ( const MathMatrix< T> &MATRIX)
	{
		// copy all elements
		if ( this != & MATRIX ) {
			// check that sizes match
			if ( NumberRows_ != MATRIX.NumberRows_ || NumberCols_ != MATRIX.NumberCols_ ) {
				// delete m_Data
				delete[] data_;

				// reallocate
				NumberRows_ = MATRIX.NumberRows_;
				NumberCols_ = MATRIX.NumberCols_;
				size_ = MATRIX.size_;
				data_ = new T[ size_ ];


			}

			// copy elements
			std::copy( MATRIX.data_, MATRIX.data_ + size_, data_);
		}
		return *this;


	}


	/// @brief assignment from value
	/// @param VALUE all elements are set to that value
	/// @return reference to this assigned Matrix
	MathMatrix< T> &operator =( const T &VALUE)
	{
		// set all element to given VALUE
		std::fill( data_, data_ + size_, VALUE);

		// return reference to this Vector
		return *this;
	}

	/// C-style data access with [] gives a pointer on a ROW
	T *operator[]( const Size ROW)
	{
		assert( ROW < NumberRows_ );
		//VectorMatrixTensorBase< t_DataType>::IsValidPosition( ROW * m_NumberCols);
		return (data_ + ROW * NumberCols_);
	}

	/// C-style data access with [] gives a pointer on a ROW
	const T *operator[]( const Size ROW) const
	{
		assert( ROW < NumberRows_ );
		//VectorMatrixTensorBase< t_DataType>::IsValidPosition( ROW * m_NumberCols);
		return (data_ + ROW * NumberCols_);
	}


	/// operator *= Matrix
	inline MathMatrix< T> &operator *= ( const MathMatrix< T> &MATRIX)
	{

		MathMatrix< T> newthis( NumberRows_, MATRIX.NumberCols_);
		for ( Size i( 0); i < newthis.NumberRows_; ++i ) {
			for ( Size j( 0); j < newthis.NumberCols_; ++j ) {
				for ( Size k( 0); k < NumberCols_; ++k ) {
					newthis( i, j) += operator()( i, k) * MATRIX( k, j);
				}
			}
		}

		return ( operator =( newthis));
	}

	/// operator *= VectorBase
	inline MathMatrix< T> & operator *= ( const MathVector< T> &VECTOR)
	{

		MathMatrix< T> newthis( NumberRows_, 1);

		const T *const vector_data( VECTOR.begin());
		for ( Size i( 0); i < NumberRows_; ++i ) {
			T sum_of_products( 0);
			T * const row = operator[]( i);
			for ( Size j( 0); j < NumberCols_; ++j ) {
				sum_of_products += row[ j] * vector_data[ j];
			}
			newthis( i, 0) = sum_of_products;
		}

		return ( operator =( newthis));
	}


private:
	Size NumberRows_; //number of rows
	Size NumberCols_; //number columns
	Size size_; //NumberRows_ * NumberCols_;
	T *data_;

};


}//end namespace numeric


#endif

