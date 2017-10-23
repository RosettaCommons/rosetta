// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/numeric.functions.hh
/// @brief  Numeric functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_numeric_functions_hh
#define INCLUDED_numeric_numeric_functions_hh


#include <numeric/types.hh>

#include <utility/assert.hh>
#include <utility/numbers.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <algorithm>
#include <cmath>
#include <limits>


namespace numeric {


// min functions: Specializations for built-in types pass by value for speed

#if (defined min) && (defined WIN32)  // Workaround for MSVC and windows.h include which used #define min
#undef min
#endif

#if (defined max) && (defined WIN32) // Workaround for MSVC and windows.h include which used #define max
#undef max
#endif

/// @brief min( short int, short int )
inline
short int
min( short int const a, short int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( int, int )
inline
int
min( int const a, int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( long int, long int )
inline
long int
min( long int const a, long int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( unsigned short int, unsigned short int )
inline
unsigned short int
min( unsigned short int const a, unsigned short int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( unsigned int, unsigned int )
inline
unsigned int
min( unsigned int const a, unsigned int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( unsigned long int, unsigned long int )
inline
unsigned long int
min( unsigned long int const a, unsigned long int const b )
{
	return ( a < b ? a : b );
}


/// @brief min( float, float )
inline
float
min( float const a, float const b )
{
	return ( a < b ? a : b );
}


/// @brief min( double, double )
inline
double
min( double const a, double const b )
{
	return ( a < b ? a : b );
}


/// @brief min( long double, long double )
inline
long double
min( long double const a, long double const b )
{
	return ( a < b ? a : b );
}


/// @brief Use std::min for 2 arguments
using std::min;


/// @brief min( a, b, c )
template< typename T >
inline
T const &
min( T const & a, T const & b, T const & c )
{
	return ( a < b ? ( a < c ? a : c ) : ( b < c ? b : c ) );
}


/// @brief min( a, b, c, d )
template< typename T >
inline
T const &
min( T const & a, T const & b, T const & c, T const & d )
{
	return std::min( std::min( a, b ), std::min( c, d ) );
}


/// @brief min( a, b, c, d, e )
template< typename T >
inline
T const &
min( T const & a, T const & b, T const & c, T const & d, T const & e )
{
	return min( std::min( a, b ), std::min( c, d ), e );
}


/// @brief min( a, b, c, d, e, f )
template< typename T >
inline
T const &
min( T const & a, T const & b, T const & c, T const & d, T const & e, T const & f )
{
	return min( std::min( a, b ), std::min( c, d ), std::min( e, f ) );
}


// max functions: Specializations for built-in types pass by value for speed


/// @brief max( short int, short int )
inline
short int
max( short int const a, short int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( int, int )
inline
int
max( int const a, int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( long int, long int )
inline
long int
max( long int const a, long int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( unsigned short int, unsigned short int )
inline
unsigned short int
max( unsigned short int const a, unsigned short int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( unsigned int, unsigned int )
inline
unsigned int
max( unsigned int const a, unsigned int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( unsigned long int, unsigned long int )
inline
unsigned long int
max( unsigned long int const a, unsigned long int const b )
{
	return ( a < b ? b : a );
}


/// @brief max( float, float )
inline
float
max( float const a, float const b )
{
	return ( a < b ? b : a );
}


/// @brief max( double, double )
inline
double
max( double const a, double const b )
{
	return ( a < b ? b : a );
}


/// @brief max( long double, long double )
inline
long double
max( long double const a, long double const b )
{
	return ( a < b ? b : a );
}


/// @brief Use std::max for 2 arguments
using std::max;


/// @brief max( a, b, c )
template< typename T >
inline
T const &
max( T const & a, T const & b, T const & c )
{
	return ( a < b ? ( b < c ? c : b ) : ( a < c ? c : a ) );
}


/// @brief max( a, b, c, d )
template< typename T >
inline
T const &
max( T const & a, T const & b, T const & c, T const & d )
{
	return std::max( std::max( a, b ), std::max( c, d ) );
}


/// @brief max( a, b, c, d, e )
template< typename T >
inline
T const &
max( T const & a, T const & b, T const & c, T const & d, T const & e )
{
	return max( std::max( a, b ), std::max( c, d ), e );
}


/// @brief max( a, b, c, d, e, f )
template< typename T >
inline
T const &
max( T const & a, T const & b, T const & c, T const & d, T const & e, T const & f )
{
	return max( std::max( a, b ), std::max( c, d ), std::max( e, f ) );
}


// General numeric functions


/// @brief square( x ) == x^2
template< typename T >
inline
T
square( T const & x )
{
	return x * x;
}


/// @brief cube( x ) == x^3
template< typename T >
inline
T
cube( T const & x )
{
	return x * x * x;
}


/// @brief sign( x )
template< typename T >
inline
int
sign( T const & x )
{
	return ( x >= T( 0 ) ? +1 : -1 );
}


/// @brief Sign transfered value
template< typename S, typename T >
inline
T
sign_transfered( S const & sigma, T const & x )
{
	return ( sigma >= S( 0 ) ? std::abs( x ) : -std::abs( x ) );
}


/// @brief Absolute difference
template< typename T >
inline
T
abs_difference( T const & a, T const & b )
{
	return max( a, b ) - min( a, b ); // Assure positive result even for unsigned types
}


/// @brief Nearest function selector class for R non-integer or T integer
template< typename R, typename T, bool >
struct NearestSelector
{
	inline
	static
	R
	nearest( T const & x )
	{
		return R( x );
	}
};


/// @brief Nearest function selector class for R integer and T non-integer
template< typename R, typename T >
struct NearestSelector< R, T, true >
{
	inline
	static
	R
	nearest( T const & x )
	{
		return R( x + ( sign( x ) * T( 0.5 ) ) );
	}
};


/// @brief nearest< R >( x ): Nearest R
template< typename R, typename T >
inline
R
nearest( T const & x )
{
	return NearestSelector< R, T, ( ( std::numeric_limits< R >::is_integer ) && ( ! std::numeric_limits< T >::is_integer ) ) >::nearest( x );
}


/// @brief nearest_size( x ): Nearest std::size_t
template< typename T >
inline
std::size_t
nearest_size( T const & x )
{
	return std::size_t( x > T( 0 ) ? x + ( sign( x ) * T( 0.5 ) ) : 0 );
}


/// @brief nearest_ssize( x ): Nearest SSize
template< typename T >
inline
SSize
nearest_ssize( T const & x )
{
	return SSize( x + ( sign( x ) * T( 0.5 ) ) );
}


/// @brief nearest_int( x ): Nearest int
template< typename T >
inline
int
nearest_int( T const & x )
{
	return static_cast< int >( x + ( sign( x ) * T( 0.5 ) ) );
}


/// @brief nint( x ): Nearest int
template< typename T >
inline
int
nint( T const & x )
{
	return static_cast< int >( x + ( sign( x ) * T( 0.5 ) ) );
}


/// @brief Mod function selector class for non-integer types
template< typename T, bool >
struct ModSelector
{
	inline
	static
	T
	mod( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( T( static_cast< SSize >( x / y ) ) * y ) : T( 0 ) );
	}
};


/// @brief Mod function selector class for integer types
/// @note When used with negative integer arguments this assumes integer division
///       rounds towards zero (de facto and future official standard)
template< typename T >
struct ModSelector< T, true >
{
	inline
	static
	T
	mod( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( ( x / y ) * y ) : T( 0 ) );
	}
};


/// @brief x(mod y) computational modulo returning magnitude < | y | and sign of x
/// @note When used with negative integer arguments this assumes integer division
///       rounds towards zero (de facto and future official standard)
template< typename T >
inline
T
mod( T const & x, T const & y )
{
	return ModSelector< T, std::numeric_limits< T >::is_integer >::mod( x, y );
}


/// @brief Modulo function selector class for non-integer types
template< typename T, bool >
struct ModuloSelector
{
	inline
	static
	T
	modulo( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( std::floor( x / y ) * y ) : T( 0 ) );
	}
};


/// @brief Modulo function selector class for integer types
template< typename T >
struct ModuloSelector< T, true >
{
	inline
	static
	T
	modulo( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( T( std::floor( static_cast< long double >( x ) / y ) ) * y ) : T( 0 ) );
	}
};


/// @brief x(mod y) mathematical modulo returning magnitude < | y | and sign of y
template< typename T >
inline
T
modulo( T const & x, T const & y )
{
	return ModuloSelector< T, std::numeric_limits< T >::is_integer >::modulo( x, y );
}


/// @brief Remainder function selector class for non-integer types
template< typename T, bool >
struct RemainderSelector
{
	inline
	static
	T
	remainder( T const & x, T const & y )
	{
		if ( y == T( 0 ) ) { // Degenerate y == 0 case
			return T( 0 );
		} else { // Normal y != 0 case
			T const x_over_y( x / y );
			SSize n( nearest_ssize( x_over_y ) );
			if ( mod( n, SSize( 2 ) ) == 1 ) { // Odd: Check for ( n - ( x / y ) ) == .5
				T const n_minus_x_over_y( T( n ) - x_over_y );
				if ( n_minus_x_over_y == T( 0.5 ) ) {
					--n;
				} else if ( n_minus_x_over_y == T( -0.5 ) ) { // Never happens if .5 rounds up
					++n;
				}
			}
			return ( x - ( n * y ) );
		}
	}
};


/// @brief Remainder function selector class for integer types
template< typename T >
struct RemainderSelector< T, true >
{
	inline
	static
	T
	remainder( T const & x, T const & y )
	{
		if ( y == T( 0 ) ) { // Degenerate y == 0 case
			return T( 0 );
		} else { // Normal y != 0 case
			long double const x_over_y( static_cast< long double >( x ) / y );
			SSize n( nearest_ssize( x_over_y ) );
			if ( mod( n, SSize( 2 ) ) == 1 ) { // Odd: Check for ( n - ( x / y ) ) == .5
				long double const n_minus_x_over_y( static_cast< long double >( n ) - x_over_y );
				if ( n_minus_x_over_y == 0.5l ) {
					--n;
				} else if ( n_minus_x_over_y == -0.5l ) { // Never happens if .5 rounds up
					++n;
				}
			}
			return ( x - ( n * y ) );
		}
	}
};


/// @brief Remainder of x with respect to division by y that is of smallest magnitude
/// @note Emulates the C99 remainder function but also supports integer arguments
/// @note Returns zero if y is zero
/// @note Return value has magnitude <= | y / 2 |
/// @note If | n - ( x / y ) | == .5 the nearest even n is used
template< typename T >
inline
T
remainder( T const & x, T const & y )
{
	return RemainderSelector< T, std::numeric_limits< T >::is_integer >::remainder( x, y );
}


/// @brief Fast remainder function selector class for non-integer types
template< typename T, bool >
struct FastRemainderSelector
{
	inline
	static
	T
	fast_remainder( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( nearest_ssize( x / y ) * y ) : T( 0 ) );
	}
};


/// @brief Fast remainder function selector class for integer types
template< typename T >
struct FastRemainderSelector< T, true >
{
	inline
	static
	T
	fast_remainder( T const & x, T const & y )
	{
		return ( y != T( 0 ) ? x - ( nearest_ssize( static_cast< long double >( x ) / y ) * y ) : T( 0 ) );
	}
};


/// @brief Remainder of x with respect to division by y that is of smallest magnitude
/// @note Emulates the C99 remainder function except for rounding halfway values to even multiples
/// @note Returns zero if y is zero
/// @note Return value has magnitude <= | y / 2 |
template< typename T >
inline
T
fast_remainder( T const & x, T const & y )
{
	return FastRemainderSelector< T, std::numeric_limits< T >::is_integer >::fast_remainder( x, y );
}


/// @brief Remainder and result of conversion to a different type
template< typename T, typename S >
inline
T
remainder_conversion( T const & t, S & s )
{
	s = S( t );
	return t - s;
}


/// @brief Greatest common divisor
template< typename T >
inline
T
gcd( T const & m, T const & n )
{
	T lo = min( m, n );
	T hi = max( m, n );
	while ( lo > T( 0 ) ) {
		T const rem = mod( hi, lo );
		hi = lo;
		lo = rem;
	}
	return hi;
}


/// @brief Equal within specified relative and absolute tolerances?
template< typename T >
inline
bool
eq_tol( T const & x, T const & y, T const & r_tol, T const & a_tol )
{
	using std::abs; // Can use std::abs or user-defined abs
	assert( r_tol >= T( 0 ) );
	assert( a_tol >= T( 0 ) );
	return ( std::abs( x - y ) <= min( r_tol * max( std::abs( x ), std::abs( y ) ), a_tol ) );
}


/// @brief Less than within specified relative and absolute tolerances?
template< typename T >
inline
bool
lt_tol( T const & x, T const & y, T const & r_tol, T const & a_tol )
{
	using std::abs; // Can use std::abs or user-defined abs
	assert( r_tol >= T( 0 ) );
	assert( a_tol >= T( 0 ) );
	return ( x < y + min( r_tol * max( std::abs( x ), std::abs( y ) ), a_tol ) );
}


/// @brief Less than or equal within specified relative and absolute tolerances?
template< typename T >
inline
bool
le_tol( T const & x, T const & y, T const & r_tol, T const & a_tol )
{
	using std::abs; // Can use std::abs or user-defined abs
	assert( r_tol >= T( 0 ) );
	assert( a_tol >= T( 0 ) );
	return ( x <= y + min( r_tol * max( std::abs( x ), std::abs( y ) ), a_tol ) );
}


/// @brief Greater than or equal within specified relative and absolute tolerances?
template< typename T >
inline
bool
ge_tol( T const & x, T const & y, T const & r_tol, T const & a_tol )
{
	using std::abs; // Can use std::abs or user-defined abs
	assert( r_tol >= T( 0 ) );
	assert( a_tol >= T( 0 ) );
	return ( x >= y - min( r_tol * max( std::abs( x ), std::abs( y ) ), a_tol ) );
}


/// @brief Greater than within specified relative and absolute tolerances?
template< typename T >
inline
bool
gt_tol( T const & x, T const & y, T const & r_tol, T const & a_tol )
{
	using std::abs; // Can use std::abs or user-defined abs
	assert( r_tol >= T( 0 ) );
	assert( a_tol >= T( 0 ) );
	return ( x > y - min( r_tol * max( std::abs( x ), std::abs( y ) ), a_tol ) );
}


/// @brief You must supply 1.0 and 0.0 as arguments a and b.
/// and no you can't short cut this, because the compiler will optimize it away!
inline
bool
is_a_finitenumber( double s, double  a, double b ){
	if ( utility::isnan(s) || utility::isinf(s) ) return false;
	if ( (a*s) != (s*cos(b)) )       return false; //  NAN!
	if ( s * 100.0 == s * 1000.00 ) return false; //  INF!
	return true;
}

// Prototype for below.
template< typename T > T factorial( T const & N );

/// @brief Calculate the value of N!.
/// @details Dangerous for large values of N.  Uses a recursive algorithm -- might not be efficient and can't be inlined.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
template< typename T >
T
factorial( T const &N ) {
	if ( N == 0 || N == 1 ) return 1;
	return N*factorial(N-1);
}


} // namespace numeric


#endif // INCLUDED_numeric_numeric_functions_HH
