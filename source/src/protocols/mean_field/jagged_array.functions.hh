// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/mean_field/jagged_array.functions.hh
/// @brief  some useful functions for jagged_array
/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)


#ifndef INCLUDED_protocols_mean_field_jagged_array_functions_HH
#define INCLUDED_protocols_mean_field_jagged_array_functions_HH


// Unit headers
#include <protocols/mean_field/jagged_array.hh>

namespace protocols {
namespace mean_field {

/// @brief Find the largest value in a jagged_array
/// @details class T must provide an operator < ()  and operator = ().
/// @details Error if input.size() == 0
template < class T >
T
max( jagged_array< T > const & input )
{
	assert( input.size() > 0 );

	T largest_so_far = utility::max ( input[ 1 ] );

	for ( typename jagged_array< T >::Size ii = 2; ii <= input.size(); ++ii ) {

		T curr_max = utility::max ( input[ ii ] );

		if ( largest_so_far < curr_max ) largest_so_far = curr_max;
	}
	return largest_so_far;
}

/// @brief Find the smallest value in a jagged_array
/// @details class T must provide an operator < ()  and operator = ().
/// @details Error if input.size() == 0
template < class T >
T
min( jagged_array< T > const & input )
{
	assert( input.size() > 0 );

	T smallest_so_far = utility::min ( input[ 1 ] );

	for ( typename jagged_array< T >::Size ii = 2; ii <= input.size(); ++ii ) {

		T curr_min = utility::min ( input[ ii ] );

		if ( smallest_so_far > curr_min ) smallest_so_far = curr_min;
	}
	return smallest_so_far;
}

/// @brief find the index of the largest value in a jagged_array
/// @details class T must provide an operator < ()  and operator = ().
/// @details Error if input.size() == 0
template < class T >
typename jagged_array< T >::Size
arg_max( jagged_array< T > const & input )
{
	assert( input.size() > 0 );

	T largest_so_far = utility::max( input[ 1 ] );

	typename jagged_array< T >::Size index_of_largest = 1;

	for ( typename jagged_array< T >::Size ii = 2; ii <= input.size(); ++ii ) {

		T curr_max = utility::max ( input[ ii ] );

		if ( largest_so_far < curr_max ) {
			largest_so_far = curr_max;
			index_of_largest = ii;
		}
	}
	return index_of_largest;
}

/// @brief find the index of the smallest value in a jagged_array
/// @details class T must provide an operator < ()  and operator = ().
/// @details Error if input.size() == 0
template < class T >
typename jagged_array< T >::Size
arg_min( jagged_array< T > const & input )
{
	assert( input.size() > 0 );

	T smallest_so_far = utility::min ( input[ 1 ] );

	typename jagged_array< T >::Size index_of_smallest = 1;

	for ( typename jagged_array< T >::Size ii = 2; ii <= input.size(); ++ii ) {

		T curr_min = utility::min ( input[ ii ] );

		if ( smallest_so_far > curr_min ) {
			smallest_so_far = curr_min;
			index_of_smallest = ii;
		}
	}
	return index_of_smallest;
}

/// @brief returns the number of elements in a jagged_array
template < class T >
inline
typename jagged_array< T >::Size
num_elements( jagged_array< T > const & input )
{
	typename jagged_array< T >::Size total = 0;

	for ( typename jagged_array< T >::Size col = 1; col <= input.size(); ++col ) {
		total += input[col].size();
	}

	return total;
}

//Operator overloading

//basic functions to be used for operator overloading

/// @brief implements add function using operator +=
/// @details class T must provide an operator +=
template < class T >
inline
T
add ( T const & operand_from_ja, T const & other_operand )
{
	T temp( operand_from_ja );
	temp += other_operand;
	return temp;
}

/// @brief implements subtract function using operator -=
/// @details class T must provide an operator -=
template < class T >
inline
T
subtract ( T const & operand_from_ja, T const & other_operand )
{
	T temp( operand_from_ja );
	temp -= other_operand;
	return temp;
}

/// @brief implements multiply function using operator *=
/// @details class T must provide an operator *=
template < class T >
inline
T
multiply ( T const & operand_from_ja, T const & other_operand )
{
	T temp( operand_from_ja );
	temp *= other_operand;
	return temp;
}

/// @brief implements divide function using operator /=
/// @details class T must provide an operator /=
template < class T >
inline
T
divide ( T const & operand_from_ja, T const & other_operand )
{
	T temp( operand_from_ja );
	temp /= other_operand;
	return temp;
}

/// @brief implements add function using operator +=
/// @details class T must provide an operator +=
/// @details allows division of class T by a different class T_2 (i.e. AAProb/Real)
template < class T, class T_2 >
inline
T
divide ( T const & operand_from_ja, T_2 const & other_operand )
{
	T temp( operand_from_ja );
	temp /= other_operand;
	return temp;
}

//Actual operator overloading
//overload three operators for each binary operation (T operand for all, utility::vector1 operands for each col, jagged_array operands for each element)
/// @brief overload operator += to allow addition of T operand to all elements in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
void
operator +=( jagged_array < T > & input, T const & operand )
{
	input.apply_func_to_all( add, operand );
}

/// @brief overload operator += to allow addition of vector1<T> operands to jagged_array in a column-dependent manner
/// @details each element in operands is added to each element of corresponding column in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
void
operator +=( jagged_array < T > & input, utility::vector1 < T > const & operands )
{
	input.apply_func_to_each_col( add, operands );
}

/// @brief overload operator += to allow addition of jagged_array<T> operands to jagged_array in an element-dependent manner
/// @details each element in operands is added to each corresponding element in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
void
operator +=( jagged_array < T > & input, jagged_array < T > const & operands )
{
	input.apply_func_to_each_elem( add, operands );
}

/// @brief overload operator + to allow addition of T operand to all elements in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
jagged_array < T >
operator +( jagged_array < T > const & input, T operand )
{
	jagged_array < T > diff (input) ;
	diff +=operand;
	return diff;
}

/// @brief overload operator + to allow addition of vector1<T> operands to jagged_array in a column-dependent manner
/// @details each element in operands is added to each element of corresponding column in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
jagged_array < T >
operator +( jagged_array < T > const & input, utility::vector1 < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff +=operands;
	return diff;
}

/// @brief overload operator + to allow addition of jagged_array<T> operands to jagged_array in an element-dependent manner
/// @details each element in operands is added to each corresponding element in jagged_array
/// @details class T must provide an operator +=
/// @details since func add uses operator +=, no need to provide operator +
template < class T >
inline
jagged_array < T >
operator +( jagged_array < T > const & input, jagged_array < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff +=operands;
	return diff;
}

/// @brief overload operator -= to allow subtraction of T operand from all elements in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
void
operator -=( jagged_array < T > & input, T operand )
{
	input.apply_func_to_all( subtract, operand );
}

/// @brief overload operator -= to allow subtraction of vector1<T> operands from jagged_array in a column-dependent manner
/// @details each element in operands is subtracted from each element of corresponding column in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
void
operator -=( jagged_array < T > & input, utility::vector1 < T > const & operands )
{
	input.apply_func_to_each_col( subtract, operands );
}


/// @brief overload operator -= to allow subtraction of jagged_array<T> operands from jagged_array in an element-dependent manner
/// @details each element in operands is subtracted from each corresponding element in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
void
operator -=( jagged_array < T > & input, jagged_array < T > const & operands )
{
	input.apply_func_to_each_elem( subtract, operands );
}

/// @brief overload operator - to allow subtraction of T operand from all elements in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
jagged_array < T >
operator -( jagged_array < T > const & input, T operand )
{
	jagged_array < T > diff (input) ;
	diff -=operand;
	return diff;
}

/// @brief overload operator - to allow subtraction of vector1<T> operands from jagged_array in a column-dependent manner
/// @details each element in operands is subtracted from each element of corresponding column in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
jagged_array < T >
operator -( jagged_array < T > const & input, utility::vector1 < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff -=operands;
	return diff;
}

/// @brief overload operator - to allow subtraction of jagged_array<T> operands from jagged_array in an element-dependent manner
/// @details each element in operands is subtracted from each corresponding element in jagged_array
/// @details class T must provide an operator -=
/// @details since func subtract uses operator -=, no need to provide operator -
template < class T >
inline
jagged_array < T >
operator -( jagged_array < T > const & input, jagged_array < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff -=operands;
	return diff;
}

/// @brief overload operator *= to allow multiplication of all elements in jagged_array by T operand
/// @details class T must provide an operator *=
/// @details since func multiply uses operator *=, no need to provide operator *
template < class T >
inline
void
operator *=( jagged_array < T > & input, T operand )
{
	input.apply_func_to_all( multiply, operand );
}

/// @brief overload operator *= to allow multiplication of jagged_array by vector1<T> operands in a column-dependent manner
/// @details each element in jagged_array is multiplied by each element of corresponding column in operands
/// @details class T must provide an operator *=
/// @details since func multiply uses operator *=, no need to provide operator *
template < class T >
inline
void
operator *=( jagged_array < T > & input, utility::vector1 < T > const & operands )
{
	input.apply_func_to_each_col( multiply, operands );
}

/// @brief overload operator *= to allow multiplication of jagged_array by jagged_array<T> operands in an element-dependent manner
/// @details each element in jagged_array is multiplied by each corresponding element in operands
/// @details class T must provide an operator *=
/// @details since func multiply uses operator *=, no need to provide operator *
template < class T >
inline
void
operator *=( jagged_array < T > & input, jagged_array < T > const & operands )
{
	input.apply_func_to_each_elem( multiply, operands );
}

/// @brief overload operator * to allow multiplication of all elements in jagged_array by T operand
/// @details class T must provide an operator *=
/// @details since func subtract uses operator *=, no need to provide operator *
template < class T >
inline
jagged_array < T >
operator *( jagged_array < T > const & input, T operand )
{
	jagged_array < T > diff (input) ;
	diff *=operand;
	return diff;
}

/// @brief overload operator * to allow multiplication of jagged_array by vector1<T> operands in a column-dependent manner
/// @details each element in jagged_array is multiplied by each element of corresponding column in operands
/// @details class T must provide an operator *=
/// @details since func multiply uses operator *=, no need to provide operator *
template < class T >
inline
jagged_array < T >
operator *( jagged_array < T > const & input, utility::vector1 < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff *=operands;
	return diff;
}

/// @brief overload operator *= to allow multiplication of jagged_array by jagged_array<T> operands in an element-dependent manner
/// @details each element in jagged_array is multiplied by each corresponding element in operands
/// @details class T must provide an operator *=
/// @details since func multiply uses operator *=, no need to provide operator *
template < class T >
inline
jagged_array < T >
operator *( jagged_array < T > const & input, jagged_array < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff *=operands;
	return diff;
}

/// @brief overload operator /= to allow division of all elements in jagged_array by T operand
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
void
operator /=( jagged_array < T > & input, T operand )
{
	input.apply_func_to_all( divide, operand );
}

/// @brief overload operator /= to allow division of jagged_array by vector1<T> operands in a column-dependent manner
/// @details each element in jagged_array is divided by each element of corresponding column in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
void
operator /=( jagged_array < T > & input, utility::vector1 < T > const & operands )
{
	input.apply_func_to_each_col( divide, operands );
}

/// @brief overload operator /= to allow division of jagged_array by vector1<T> operands in a column-dependent manner
/// @details each element in jagged_array is divided by each element of corresponding column in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
/// @details intended to allow for division of a class T by a different class T_2 - used in mean-field method to divide AAProb/Real
template < class T, class T_2 >
inline
void
operator /=( jagged_array < T > & input, utility::vector1 < T_2 > const & operands )
{
	input.apply_func_to_each_col( divide, operands );
}

/// @brief overload operator /= to allow division of jagged_array by jagged_array<T> operands in an element-dependent manner
/// @details each element in jagged_array is divided by each each corresponding element in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
void
operator /=( jagged_array < T > & input, jagged_array < T > const & operands )
{
	input.apply_func_to_each_elem( divide, operands );
}

/// @brief overload operator /= to allow division of jagged_array by jagged_array<T> operands in an element-dependent manner
/// @details each element in jagged_array is divided by each each corresponding element in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
/// @details intended to allow for division of a class T by a different class T_2 - used in mean-field method to divide AAProb/Real
template < class T, class T_2 >
inline
void
operator /=( jagged_array < T > & input, jagged_array < T_2 > const & operands )
{
	input.apply_func_to_each_elem( divide, operands );
}

/// @brief overload operator / to allow division of all elements in jagged_array by T operand
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
jagged_array < T >
operator /( jagged_array < T > const & input, T operand )
{
	jagged_array < T > diff (input) ;
	diff /=operand;
	return diff;
}

/// @brief overload operator / to allow division of jagged_array by vector1<T> operands in a column-dependent manner
/// @details each element in jagged_array is divided by each element of corresponding column in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
jagged_array < T >
operator /( jagged_array < T > const & input, utility::vector1 < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff /=operands;
	return diff;
}

/// @brief overload operator / to allow division of jagged_array by jagged_array<T> operands in an element-dependent manner
/// @details each element in jagged_array is divided by each each corresponding element in operands
/// @details class T must provide an operator /=
/// @details since func divide uses operator /=, no need to provide operator /
template < class T >
inline
jagged_array < T >
operator /( jagged_array < T > const & input, jagged_array < T > const & operands )
{
	jagged_array < T > diff (input) ;
	diff /=operands;
	return diff;
}

} // namespace mean_field
} // namespace protocols

#endif
