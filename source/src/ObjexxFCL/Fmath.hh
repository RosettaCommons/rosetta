#ifndef INCLUDED_ObjexxFCL_Fmath_hh
#define INCLUDED_ObjexxFCL_Fmath_hh


// Fortran Intrinsic-Compatible and General Math Functions
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// C++ Headers
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>


namespace ObjexxFCL {


typedef  ptrdiff_t  SSize;


// min Functions


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


// max Functions


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


// Math Functions


/// @brief std::abs( x ) == | x |
template< typename T >
inline
T
abs( T const & x )
{
	return ( x < T( 0 ) ? -x : x );
}


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


/// @brief Sign Transfer from Second Argument to First Argument
template< typename X, typename Y >
inline
X
sign( X const & x, Y const & y )
{
	return ( y >= Y( 0 ) ? std::abs( x ) : -std::abs( x ) );
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


/// @brief nsint( x ): Nearest short int
template< typename T >
inline
short int
nsint( T const & x )
{
	return static_cast< short int >( x + ( sign( x ) * T( 0.5 ) ) );
}


/// @brief nlint( x ): Nearest long int
template< typename T >
inline
long int
nlint( T const & x )
{
	return static_cast< long int >( x + ( sign( x ) * T( 0.5 ) ) );
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


/// @brief i(mod n) : float Arguments
inline
float
mod( float const & i, float const & n )
{
	assert( n != 0.0f );
	return ( n == 0.0f ? 0.0f : std::fmod( i, n ) );
}


/// @brief i(mod n) : double Arguments
inline
double
mod( double const & i, double const & n )
{
	assert( n != 0.0 );
	return ( n == 0.0 ? 0.0 : std::fmod( i, n ) );
}


/// @brief i(mod n) : long double Arguments
inline
long double
mod( long double const & i, long double const & n )
{
	assert( n != 0.0l );
	return ( n == 0.0l ? 0.0l : std::fmod( i, n ) );
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


/// @brief Greatest Common Divisor
template< typename T >
inline
T
gcd( T const & m, T const & n )
{
	T lo( min( m, n ) );
	T hi( max( m, n ) );
	while ( lo > T( 0 ) ) {
		T const rem( mod( hi, lo ) );
		hi = lo;
		lo = rem;
	}
	return hi;
}


// Comparison-With-Tolerance Functions


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


// Bit Functions
namespace bit { // Protect from collisions with C++0x <functional> functors


/// @brief Bitwise Not
template< typename T >
inline
T
bit_not( T const & x )
{
	return ( ~x );
}


/// @brief Bitwise And
template< typename T >
inline
T
bit_and( T const & x, T const & y )
{
	return ( x & y );
}


/// @brief Bitwise Inclusive Or
template< typename T >
inline
T
bit_or( T const & x, T const & y )
{
	return ( x | y );
}


/// @brief Bitwise Exclusive Or
template< typename T >
inline
T
bit_xor( T const & x, T const & y )
{
	return ( x ^ y );
}


/// @brief Bit Value Set to 1
template< typename T >
inline
T
bit_set( T const & x, T const & pos )
{
	assert( pos >= T( 0 ) );
	return ( x | ( T( 1 ) << pos ) );
}


/// @brief Bit Value Set to 0
template< typename T >
inline
T
bit_clr( T const & x, T const & pos )
{
	assert( pos >= T( 0 ) );
	return ( x & ~( T( 1 ) << pos ) );
}


/// @brief Bit Value Test
template< typename T >
inline
bool
bit_test( T const & x, T const & pos )
{
	assert( pos >= T( 0 ) );
	return ( x & ( T( 1 ) << pos ) );
}


} // namespace bit


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_Fmath_HH
