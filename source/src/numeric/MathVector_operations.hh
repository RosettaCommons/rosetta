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
/// Mathmatical functions for the MathMatrix class
///
/// @details
/// This is an implementation of an algorithm that was taken from BCL (Jens Meiler)
/// The Matrix is construction is found in numeric/MathVector.hh. These are mathematical functions
/// that can be used by the MathVector class. ***Note that these are outside of class, but having the
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


#ifndef INCLUDED_numeric_MathVector_operations_hh
#define INCLUDED_numeric_MathVector_operations_hh

#include <numeric/MathVector.hh>

namespace numeric{


//////////////////////
// Vector functions //
//////////////////////


/// @return distance between two MathVectors A->B
template< typename T >
inline T distance( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return( ( VECTOR_B - VECTOR_A).norm());
}

/// @return  projection angle between four  MathVectors A->B and C->D
template< typename T >
inline T proj_angl
(
	MathVector< T > const & VECTOR_A,
	MathVector< T > const & VECTOR_B,
	MathVector< T > const & VECTOR_C,
	MathVector< T > const & VECTOR_D
)
{
	MathVector< T > ab = VECTOR_B - VECTOR_A, cd = VECTOR_D - VECTOR_C;
	T angle = (ab * cd) / ( ab.norm() * cd.norm());
	angle = acos( std::min( T( 1), std::max( T( -1), angle)));
	return angle;
}

/// @return projection angle between three MathVectors A->B and A->C
template< typename T >
inline T proj_angl (
	MathVector< T > const & VECTOR_A,
	MathVector< T > const & VECTOR_B,
	MathVector< T > const & VECTOR_C
)
{
	return proj_angl( VECTOR_A, VECTOR_B, VECTOR_A, VECTOR_C);
}

/// @return projection angle between two   MathVectors 0->A and 0->B
template< typename T >
inline T proj_angl( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return proj_angl( MathVector< T >( 3), VECTOR_A, MathVector< T >( 3), VECTOR_B);
}

/// @return scalar product of two MathVectors
template< typename T >
T scalar_product( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return ( VECTOR_A * VECTOR_B);
}

template< typename T >
MathVector< T > MakeVector( T const & X)
{
	return MathVector< T >( 1, X);
}

/// @return construct MathVectors from two elements
template< typename T >
inline MathVector< T > MakeVector( T const & X, T const & Y)
{
	MathVector< T > newvector( 2);
	newvector( 0) = X;
	newvector( 1) = Y;

	return newvector;
}

/// @return construct MathVectors from three elements
template< typename T >
inline	MathVector< T > MakeVector( T const & X, T const & Y, T const & Z)
{
	MathVector< T > newvector( 3 );
	newvector( 0 ) = X;
	newvector( 1 ) = Y;
	newvector( 2 ) = Z;

	return newvector;
}


///////////////
// operators //
///////////////

/// @return operator -MathVectors
template< typename T >
inline MathVector< T > operator - ( MathVector< T > const & VECTOR)
{
	return MathVector< T >( T( -1) * VECTOR);
}

/// @return operator +MathVectors
template< typename T >
inline MathVector< T > operator + ( MathVector< T > const & VECTOR)
{
	return MathVector< T >( VECTOR);
}

/// @return operator == (Comparison) MathVectors
template< typename T >
inline bool operator ==( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	if( VECTOR_A.size() != VECTOR_B.size()) {
		return false;
	}
	for ( const T
			*ptr_a( VECTOR_A.begin()), *ptr_a_end( VECTOR_A.end()),
			*ptr_b = VECTOR_B.begin(), *ptr_b_end( VECTOR_B.end());
			ptr_a != ptr_a_end && ptr_b != ptr_b_end; ++ptr_a, ++ptr_b ) {
		if( ( *ptr_a) != ( *ptr_b)) {
			return false;
		}
	}
	return true;
}

/// @return operator != (Comparison) MathVectors
template< typename T >
inline bool operator != ( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return !( VECTOR_A == VECTOR_B);
}

/// @return operator MathVector == X (Comparison with T value) MathVectors
template< typename T >
inline bool operator ==( MathVector< T > const & VECTOR, T const & X)
{
	for( const T *ptr( VECTOR.begin()), *ptr_end( VECTOR.end()); ptr != ptr_end; ++ptr) {
		if( ( *ptr) != X) return false;
	}
	return true;
}

/// @return operator MathVector != X (Comparison with T value)
template< typename T >
inline bool operator !=( MathVector< T > const & VECTOR, T const & X)
{
	return !( VECTOR == X);
}

/// @return operator X == MathVector(Comparison with T value)
template< typename T >
inline bool operator ==( T const & X, MathVector< T > const & VECTOR)
{
	return( VECTOR == X);
}

/// @return operator X != MathVector(Comparison with T value)
template< typename T >
inline bool operator !=( T const & X, MathVector< T > const & VECTOR)
{
	return !( VECTOR == X);
}

/// @return operator MathVector + MathVector
template< typename T >
inline MathVector< T > operator + ( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return MathVector< T >( VECTOR_A).operator +=( VECTOR_B);
}

/// @return operator MathVector - MathVector
template< typename T >
inline MathVector< T > operator - ( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	return MathVector< T >( VECTOR_A ).operator -=( VECTOR_B);
}

/// @return operator MathVector + MathVector
template< typename T >
inline MathVector< T > operator + ( MathVector< T > const & VECTOR, T const & X)
{
	return MathVector< T >( VECTOR ).operator +=( X);
}

/// @return  operator MathVector - MathVector
template< typename T >
inline MathVector< T > operator - ( MathVector< T > const & VECTOR, T const & X)
{
	return MathVector< T >( VECTOR ).operator -=( X);
}

/// @return operator MathVector + MathVector
template< typename T >
inline MathVector< T > operator + ( T const & X, MathVector< T > const & VECTOR)
{
	return MathVector< T >( VECTOR ).operator +=( X);
}

/// @return operator MathVector - MathVector
template< typename T >
inline MathVector< T > operator - ( T const & X, MathVector< T > const & VECTOR)
{
	return MathVector< T >( -VECTOR ).operator +=( X);
}

/// @return operator MathVector * MathVector
template< typename T >
inline T operator * ( MathVector< T > const & VECTOR_A, MathVector< T > const & VECTOR_B)
{
	T scalar( 0);
	for ( const T
			*ptr_a( VECTOR_A.begin()), *ptr_a_end( VECTOR_A.end()),
			*ptr_b( VECTOR_B.begin()), *ptr_b_end( VECTOR_B.end());
			ptr_a != ptr_a_end && ptr_b != ptr_b_end; ++ptr_a, ++ptr_b ) {
		scalar += ( *ptr_a) * ( *ptr_b);
	}

	return scalar;
}

/// @return operator MathVector * MathVector
template< typename T >
inline MathVector< T > operator *( T const & X, MathVector< T > const & VECTOR)
{
	return MathVector< T >( VECTOR).operator *=( X);
}

/// @return operator MathVector * MathVector
template< typename T >
inline MathVector< T > operator *( MathVector< T > const & VECTOR,  T const & X)
{
	return MathVector< T >( VECTOR).operator *=( X);
}

/// @return operator MathVector / MathVector
template< typename T >
inline MathVector< T > operator /( MathVector< T > const & VECTOR,  T const & X)
{
	return MathVector< T >( VECTOR).operator /=( X);
}

/// @return operator MathVector / MathVector ; divides each element by according argument element
template< typename T >
inline MathVector< T > operator /( MathVector< T > const & VECTOR_A,  MathVector< T > const & VECTOR_B)
{
	return MathVector< T >( VECTOR_A).operator /=( VECTOR_B);
}

/// @return operator MathVector / MathVector = Value * Inverse( MathVector)
template< typename T >
inline MathVector< T > operator / ( T const & X, MathVector< T > const & VECTOR)
{
	return MathVector< T >( Inverse( VECTOR)).operator *=( X);
}

/// @return MathVector operator Value ^ MathVector
template< typename T >
inline MathVector< T > operator ^ ( T const & X, MathVector< T > const & VECTOR)
{
	MathVector< T > newvector( VECTOR.GetSize());
	T *ptr_new( newvector.begin()), *ptr_new_end( newvector.end());
	for ( const T
			*ptr_old( VECTOR.begin()), *ptr_old_end( VECTOR.end());
			ptr_old != ptr_old_end && ptr_new != ptr_new_end;
			++ptr_new, ++ptr_old ) {
		( *ptr_new) = pow( X, ( *ptr_old));
	}
	return newvector;
}


}

#endif /* MATHVECTOR_OPERATIONS_HH_ */
