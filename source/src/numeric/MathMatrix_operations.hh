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
/// Mathmatical functions for the MathMatrix class
///
/// @details
/// This is an implementation of an algorithm that was taken from BCL (Jens Meiler)
/// The Matrix is construction is found in numeric/MathMatrix.hh. These are mathematical functions
/// that can be used by the MathMatrix class. ***Note that these are outside of class, but having the
/// operators take two arguments ensures that no one will use the functions unknowningly
///
///
/// @references
/// Nils Woetzl
/// Jens Meiler
///
/// @author Steven Combs, Nils Woetzl, Jens Meiler
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathMatrix_operations_hh
#define INCLUDED_numeric_MathMatrix_operations_hh

#include <numeric/MathMatrix.hh>

#include <algorithm> //needed for std::transform, std::find_if
namespace numeric {


//////////////////////
// binary operators //
//////////////////////

/// @brief add one matrix to another
/// @param MATRIX_LHS matrix to add to
/// @param MATRIX_RHS matrix to add
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator += ( MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{

	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_RHS.begin(),
		MATRIX_LHS.begin(),
		std::plus< T>()
	);

	return MATRIX_LHS;
}

/// @brief subtract one matrix from another
/// @param MATRIX_LHS matrix to subtract from
/// @param MATRIX_RHS matrix to subtract
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator -= ( MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{


	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_RHS.begin(),
		MATRIX_LHS.begin(),
		std::minus< T>()
	);

	return MATRIX_LHS;
}

/// @brief divide one matrix by another
/// @param MATRIX_LHS matrix to divided
/// @param MATRIX_RHS matrix to divide by
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator /= ( MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> & )
{

	return MATRIX_LHS;
}

/// @brief add scalar to matrix
/// @param MATRIX_LHS matrix to add to
/// @param VALUE scalar to be added
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator += ( MathMatrix< T> &MATRIX_LHS, const T &VALUE )
{
	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_LHS.begin(),
		std::bind2nd( std::plus< T>(), VALUE)
	);

	return MATRIX_LHS;
}

/// @brief subtract scalar from matrix
/// @param MATRIX_LHS matrix to subtract from
/// @param VALUE scalar to be added
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator -= ( MathMatrix< T> &MATRIX_LHS, const T &VALUE )
{
	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_LHS.begin(),
		std::bind2nd( std::minus< T>(), VALUE)
	);

	return MATRIX_LHS;
}

/// @brief multiply matrix with scalar
/// @param MATRIX_LHS matrix to multiply to
/// @param SCALAR scalar to be multiplied
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator *= ( MathMatrix< T> &MATRIX_LHS, const T &SCALAR )
{
	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_LHS.begin(),
		std::bind2nd( std::multiplies< T>(), SCALAR)
	);

	return MATRIX_LHS;
}

/// @brief divide matrix by scalar
/// @param MATRIX_LHS matrix to divide
/// @param SCALAR scalar to divide by
/// @return the changed lhs matrix
template< typename T>
inline
MathMatrix< T> & operator /= ( MathMatrix< T> &MATRIX_LHS, const T &SCALAR )
{
	std::transform
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		MATRIX_LHS.begin(),
		std::bind2nd( std::divides< T>(), SCALAR)
	);

	return MATRIX_LHS;
}

//////////////////////////////
// binary logical operators //
//////////////////////////////

/// @brief compare to matricess for equality
/// @param MATRIX_LHS lhs matrix
/// @param MATRIX_RHS rhs matrix
/// @return true is they are equal in size and all pairs of items are equal
template< typename T>
inline
bool operator == ( const MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	// not equal if different size
	if ( !same_dimensions( MATRIX_LHS, MATRIX_RHS) ) {
		return false;
	}

	// check that all elements in bot ranges comparing each corresponding pair, is equal
	return std::equal( MATRIX_LHS.begin(), MATRIX_LHS.end(), MATRIX_RHS.begin());
}

/// @brief compare to matrices for inequality
/// @param MATRIX_LHS lhs matrix
/// @param MATRIX_RHS rhs matrix
/// @return !( MATRIX_LHS == MATRIX_RHS)
template< typename T>
inline
bool operator != ( const MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	return !( MATRIX_LHS == MATRIX_RHS);
}

/// @brief compare if all items in matrix are equal to a given VALUE
/// @param MATRIX_LHS matrix with values
/// @param VALUE_RHS value that is compared against
/// @return true if matrix is empty are all elements in matrix are equal to given VALUE
template< typename T>
inline
bool operator == ( const MathMatrix< T> &MATRIX_LHS, const T &VALUE_RHS )
{
	return std::find_if
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		std::bind2nd( std::not_equal_to< T>(), VALUE_RHS)
		) == MATRIX_LHS.end();
}

/// @brief compare if all items in matrix are equal to a given VALUE
/// @param VALUE_LHS value that is compared against
/// @param MATRIX_RHS matrix with values
/// @return true if matrix is empty are all elements in matrix are equal to given VALUE
template< typename T>
inline
bool operator == ( const T &VALUE_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	return ( MATRIX_RHS == VALUE_LHS);
}

/// @brief compare if all items in matrix are not equal to a given VALUE
/// @param MATRIX_LHS matrix with values
/// @param VALUE_RHS value that is compared against
/// @return false if matrix is empty are all elements in matrix are equal to given VALUE
template< typename T>
inline
bool operator != ( const MathMatrix< T> &MATRIX_LHS, const T &VALUE_RHS )
{
	return std::find_if
		(
		MATRIX_LHS.begin(), MATRIX_LHS.end(),
		std::bind2nd( std::not_equal_to< T>(), VALUE_RHS)
		) != MATRIX_LHS.end();
}

/// @brief compare if all items in matrix are not equal to a given VALUE
/// @param VALUE_LHS value that is compared against
/// @param MATRIX_RHS matrix with values
/// @return false if matrix is empty are all elements in matrix are equal to given VALUE
template< typename T>
inline
bool operator != ( const T &VALUE_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	return ( MATRIX_RHS != VALUE_LHS);
}

//////////////////////
// binary operators //
//////////////////////

/// @brief sum two matrixs of equal size
/// @param MATRIX_LHS lhs matrix
/// @param MATRIX_RHS rhs matrix
/// @return matrix with all individual summed elements of lhs and rhs matrix
template< typename T>
inline
MathMatrix< T>
operator + ( const MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{


	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix += MATRIX_RHS;
}

/// @brief subtract two matrixs of equal size
/// @param MATRIX_LHS lhs matrix
/// @param MATRIX_RHS rhs matrix
/// @return matrix with all individual subtracted elements of rhs from lhs matrix
template< typename T>
inline
MathMatrix< T> operator - ( const MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{


	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix -= MATRIX_RHS;
}


/// @brief multiply two matrixs of equal size by building the inner product yielding the scalar product
/// @param MATRIX_LHS lhs matrix
/// @param MATRIX_RHS rhs matrix
/// @return scalar representing root of inner product of the two ranges
template< typename T>
inline
MathMatrix< T> operator * ( const MathMatrix< T> &MATRIX_LHS, const MathMatrix< T> &MATRIX_RHS )
{


	MathMatrix< T> new_matrix( MATRIX_LHS.get_number_rows(), MATRIX_RHS.get_number_cols());
	for ( Size i( 0); i < new_matrix.get_number_rows(); ++i ) {
		for ( Size j( 0); j < new_matrix.get_number_cols(); ++j ) {
			for ( Size k( 0); k < MATRIX_LHS.get_number_cols(); ++k ) {
				new_matrix( i, j) += MATRIX_LHS.operator()( i, k) * MATRIX_RHS( k, j);
			}
		}
	}
	return new_matrix;
}

/// @brief add value to matrix
/// @param MATRIX_LHS lhs matrix
/// @param VALUE_RHS rhs value to be added
/// @return matrix that has the value added to each value of the lhs given matrix
template< typename T>
inline
MathMatrix< T> operator + ( const MathMatrix< T> &MATRIX_LHS, const T &VALUE_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix += VALUE_RHS;
}

/// @brief add matrix to value
/// @param VALUE_LHS lhs value to be added
/// @param MATRIX_RHS rhs matrix
/// @return matrix that has the value added to each value of the lhs given matrix
template< typename T>
inline
MathMatrix< T> operator + ( const T &VALUE_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_RHS);
	return new_matrix += VALUE_LHS;
}

/// @brief subtract value from matrix
/// @param MATRIX_LHS lhs matrix
/// @param VALUE_RHS rhs value to be subtracted
/// @return matrix that has the value subtracted from each value of the lhs given matrix
template< typename T>
inline
MathMatrix< T> operator - ( const MathMatrix< T> &MATRIX_LHS, const T &VALUE_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix -= VALUE_RHS;
}

/// @brief subtract matrix from value
/// @param VALUE_LHS rhs value to be subtracted
/// @param MATRIX_RHS lhs matrix
/// @return matrix that has the values in the matrix subtracted from the value
template< typename T>
inline
MathMatrix< T> operator - ( const T &VALUE_LHS, const MathMatrix< T> &MATRIX_RHS)
{
	MathMatrix< T> new_matrix( MATRIX_RHS.size(), VALUE_LHS);
	return new_matrix -= MATRIX_RHS;
}

/// @brief multiply scalar with matrix
/// @param SCALAR_LHS lhs value to be multiplied
/// @param MATRIX_RHS rhs matrix
/// @return matrix that has the values multiplied with the scalar
template< typename T>
inline
MathMatrix< T> operator * ( const T &SCALAR_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_RHS);
	return new_matrix *= SCALAR_LHS;
}

/// @brief multiply matrix with scalar
/// @param MATRIX_LHS lhs matrix
/// @param SCALAR_RHS rhs value to be multiplied
/// @return matrix that has the values multiplied with the scalar
template< typename T>
inline
MathMatrix< T> operator * ( const MathMatrix< T> &MATRIX_LHS, const T &SCALAR_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix *= SCALAR_RHS;
}

/// @brief multiply matrix with vector
/// @param MATRIX_LHS lhs matrix
/// @param VECTOR vector to be multiplied
/// @return resulting vector
template< typename T>
inline
MathVector< T> operator * ( const MathMatrix< T> &MATRIX_LHS, const MathVector< T> &VECTOR_RHS )
{
	MathVector< T> new_vector( MATRIX_LHS.get_number_rows());
	for ( Size i( 0); i < MATRIX_LHS.get_number_rows(); ++i ) {
		for ( Size j( 0); j < MATRIX_LHS.get_number_cols(); ++j ) {
			new_vector( i ) += MATRIX_LHS.operator()( i, j) * VECTOR_RHS( j);
		}
	}

	return new_vector;
}

/// @brief divide matrix with scalar
/// @param MATRIX_LHS lhs matrix
/// @param SCALAR_RHS rhs value to be divided by
/// @return matrix that has the values divided by the scalar
template< typename T>
inline
MathMatrix< T> operator / ( const MathMatrix< T> &MATRIX_LHS, const T &SCALAR_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_LHS);
	return new_matrix /= SCALAR_RHS;
}

/// @brief divide scalar by matrix
/// @param SCALAR_LHS lhs value to be divided
/// @param MATRIX_RHS rhs matrix to be used to divide the scalar
/// @return matrix that has the values of scalar divided by each according value of the matrix
template< typename T>
inline
MathMatrix< T> operator / ( const T &SCALAR_LHS, const MathMatrix< T> &MATRIX_RHS )
{
	MathMatrix< T> new_matrix( MATRIX_RHS.size(), SCALAR_LHS);
	return new_matrix /= MATRIX_RHS;
}

}


#endif /* MATHMATRIX_OPERATIONS_HH_ */


