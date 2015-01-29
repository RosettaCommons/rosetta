// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/utility.functions.hh
/// @brief  Numeric functions (to avoid dependency on numeric package)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_utility_functions_hh
#define INCLUDED_utility_utility_functions_hh


// C++ headers
#include <algorithm>
#include <utility/assert.hh>
#include <cmath>
#include <limits>


namespace utility {


// min functions: Specializations for built-in types pass by value for speed


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


} // namespace utility


#endif // INCLUDED_utility_utility_functions_HH
