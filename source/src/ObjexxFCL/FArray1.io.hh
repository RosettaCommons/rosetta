#ifndef INCLUDED_ObjexxFCL_FArray1_io_hh
#define INCLUDED_ObjexxFCL_FArray1_io_hh


// FArray1 Input/Output Functions
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


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/TypeTraits.hh>

// C++ Headers
#include <cstddef>
#include <iomanip>
#include <istream>
#include <ostream>


namespace ObjexxFCL {


/// @brief Read an FArray1 from a binary file
template< typename T >
std::istream &
read_binary( std::istream & stream, FArray1< T > & a )
{
	// Read array from stream in column-major (Fortran) order
	if ( stream ) {
		std::size_t const type_size( sizeof( T ) / sizeof( std::istream::char_type ) );
		for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
			if ( stream ) stream.read( ( std::istream::char_type * )&a( i ), type_size );
		}
	}
	return stream;
}


/// @brief Write an FArray1 to a binary file
template< typename T >
std::ostream &
write_binary( std::ostream & stream, FArray1< T > const & a )
{
	// Write array to stream in column-major (Fortran) order
	if ( stream ) {
		std::size_t const type_size( sizeof( T ) / sizeof( std::ostream::char_type ) );
		for ( int i = a.l(), e = a.u(); i <= e; ++i ) {
			if ( stream ) stream.write( ( std::ostream::char_type const * )&a( i ), type_size );
		}
	}
	return stream;
}


/// @brief stream << FArray1
template< typename T >
std::ostream &
operator <<( std::ostream & stream, FArray1< T > const & a )
{
	if ( a.size() == 0 ) return stream;

	if ( stream ) {

		// Types
		using std::setw;
		typedef  TypeTraits< T >  Traits;

		// Save current stream state and set persistent state
		std::ios_base::fmtflags const old_flags( stream.flags() );
		int const old_precision( stream.precision( Traits::precision() ) );
		stream << std::right << std::showpoint << std::uppercase;

		// Output array to stream
		int const w( Traits::width() );
		for ( int i = a.l(), e = a.u(); i < e; ++i ) {
			stream << setw( w ) << a( i ) << ' ';
		} stream << setw( w ) << a( a.u() );

		// Restore previous stream state
		stream.precision( old_precision );
		stream.flags( old_flags );

	}

	return stream;
}


namespace format {


/// @brief Standard Width Format: FArray1
template< typename T >
inline
std::string
SW( FArray1< T > const & a )
{
	std::size_t const n( a.size() );
	std::string s;
	s.reserve( n * static_cast< std::size_t >( TypeTraits< T >::standard_width() ) );
	for ( std::size_t i = 0; i < n; ++i ) {
		s.append( format::SW( a[ i ] ) );
	}
	return s;
}


} // namespace format


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1_io_HH
