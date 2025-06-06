#ifndef INCLUDED_ObjexxFCL_FArray6_io_hh
#define INCLUDED_ObjexxFCL_FArray6_io_hh


// FArray6 Input/Output Functions
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
#include <ObjexxFCL/FArray6.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>

// C++ Headers
#include <iomanip>
#include <istream>


namespace ObjexxFCL {


/// @brief Read an FArray6 from a binary file
template< typename T >
std::istream &
read_binary( std::istream & stream, FArray6< T > & a )
{
	// Read array from stream in column-major (Fortran) order
	if ( stream ) {
		std::size_t const type_size( sizeof( T ) / sizeof( std::istream::char_type ) );
		for ( int i6 = a.l6(), e6 = a.u6(); i6 <= e6; ++i6 ) {
			if ( stream ) {
				for ( int i5 = a.l5(), e5 = a.u5(); i5 <= e5; ++i5 ) {
					if ( stream ) {
						for ( int i4 = a.l4(), e4 = a.u4(); i4 <= e4; ++i4 ) {
							if ( stream ) {
								for ( int i3 = a.l3(), e3 = a.u3(); i3 <= e3; ++i3 ) {
									if ( stream ) {
										for ( int i2 = a.l2(), e2 = a.u2(); i2 <= e2; ++i2 ) {
											if ( stream ) {
												for ( int i1 = a.l1(), e1 = a.u1(); i1 <= e1; ++i1 ) {
													if ( stream ) stream.read( ( std::istream::char_type * )&a( i1, i2, i3, i4, i5, i6 ), type_size );
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return stream;
}


/// @brief Write an FArray6 to a binary file
template< typename T >
std::ostream &
write_binary( std::ostream & stream, FArray6< T > const & a )
{
	// Write array to stream in column-major (Fortran) order
	if ( stream ) {
		std::size_t const type_size( sizeof( T ) / sizeof( std::ostream::char_type ) );
		for ( int i6 = a.l6(), e6 = a.u6(); i6 <= e6; ++i6 ) {
			if ( stream ) {
				for ( int i5 = a.l5(), e5 = a.u5(); i5 <= e5; ++i5 ) {
					if ( stream ) {
						for ( int i4 = a.l4(), e4 = a.u4(); i4 <= e4; ++i4 ) {
							if ( stream ) {
								for ( int i3 = a.l3(), e3 = a.u3(); i3 <= e3; ++i3 ) {
									if ( stream ) {
										for ( int i2 = a.l2(), e2 = a.u2(); i2 <= e2; ++i2 ) {
											if ( stream ) {
												for ( int i1 = a.l1(), e1 = a.u1(); i1 <= e1; ++i1 ) {
													if ( stream ) stream.write( ( std::ostream::char_type const * )&a( i1, i2, i3, i4, i5, i6 ), type_size );
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return stream;
}


/// @brief stream << FArray6
template< typename T >
std::ostream &
operator <<( std::ostream & stream, FArray6< T > const & a )
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
		for ( int i1 = a.l1(), e1 = a.u1(); i1 <= e1; ++i1 ) {
			for ( int i2 = a.l2(), e2 = a.u2(); i2 <= e2; ++i2 ) {
				for ( int i3 = a.l3(), e3 = a.u3(); i3 <= e3; ++i3 ) {
					for ( int i4 = a.l4(), e4 = a.u4(); i4 <= e4; ++i4 ) {
						for ( int i5 = a.l5(), e5 = a.u5(); i5 < e5; ++i5 ) {
							for ( int i6 = a.l6(), e6 = a.u6(); i6 < e6; ++i6 ) {
								stream << setw( w ) << a( i1, i2, i3, i4, i5, i6 ) << ' ';
							} stream << setw( w ) << a( i1, i2, i3, i4, i5, a.u6() ) << '\n';
						}
					}
				}
			}
		}

		// Restore previous stream state
		stream.precision( old_precision );
		stream.flags( old_flags );

	}

	return stream;
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray6_io_HH
