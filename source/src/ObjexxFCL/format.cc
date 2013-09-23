// Fortran-Compatible Formatted Input/Output Support
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
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/byte.hh>
#include <ObjexxFCL/ubyte.hh>
#include <ObjexxFCL/Fstring.hh>

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <iostream>


namespace ObjexxFCL {
namespace format {


// Formatted Input


// Bite Explicit Specializations


	/// @brief Assign Stream Bite to Value: byte Specialization
	template<>
	void
	Bite< byte >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: ubyte Specialization
	template<>
	void
	Bite< ubyte >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: Fstring Specialization
	template<>
	void
	Bite< Fstring >::assign( std::stringstream & ss ) const
	{
		t_ = ss.str();
	}


// Bite Makers


/// @brief Fstring Bite Maker: Take Length of Fstring
Bite< Fstring >
bite( Fstring & t )
{
	return Bite< Fstring >( t.length(), t );
}


// Formatted Output


// General Formatting


/// @brief Single-Spaced Format: bool Specialization
template<>
std::string
SS( bool const & t )
{
	return L( 2, t );
}


/// @brief Single-Spaced Format: float Specialization
template<>
std::string
SS( float const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << ( t < 0.0f ? " " : "  " ) << t;
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: double Specialization
template<>
std::string
SS( double const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << ( t < 0.0 ? " " : "  " ) << t;
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: long double Specialization
template<>
std::string
SS( long double const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << ( t < 0.0l ? " " : "  " ) << t;
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: complex< float > Specialization
template<>
std::string
SS( std::complex< float > const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << " ("
	 << ( t.real() < 0.0f ? "" : " " ) << t.real() << ','
	 << ( t.imag() < 0.0f ? "" : " " ) << t.imag() << ')';
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: complex< double > Specialization
template<>
std::string
SS( std::complex< double > const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << " ("
	 << ( t.real() < 0.0 ? "" : " " ) << t.real() << ','
	 << ( t.imag() < 0.0 ? "" : " " ) << t.imag() << ')';
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: complex< long double > Specialization
template<>
std::string
SS( std::complex< long double > const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << std::setprecision( 9 ) << " ("
	 << ( t.real() < 0.0l ? "" : " " ) << t.real() << ','
	 << ( t.imag() < 0.0l ? "" : " " ) << t.imag() << ')';
	return fmt_stream.str();
}


// String Formatting


/// @brief char Format
std::string
A( int const w, char const c )
{
	if ( w <= 0 ) {
		return std::string();
	} else if ( w == 1 ) {
		return std::string( 1, c );
	} else {
		return std::string( w-1, ' ' ) + c;
	}
}


/// @brief char Format (Width==1)
std::string
A( char const c )
{
	return std::string( 1, c );
}


/// @brief cstring Format
std::string
A( int const w, c_cstring const s )
{
	return A( w, std::string( s ) );
}


/// @brief cstring Format (Width of cstring)
std::string
A( c_cstring const s )
{
	return std::string( s );
}


/// @brief string Format
std::string
A( int const w, std::string const & s )
{
	std::string::size_type const s_length( s.length() );
	if ( w <= 0 ) {
		return std::string();
	} else if ( static_cast< int >( s_length ) > w ) {
		return s.substr( 0, w );
	} else if ( static_cast< int >( s_length ) == w ) {
		return s;
	} else { // s_length < w
		return std::string( w - s_length, ' ' ) + s;
	}
}


/// @brief string Format (Width of string)
std::string const &
A( std::string const & s )
{
	return s;
}


/// @brief Fstring Format
std::string
A( int const w, Fstring const & s )
{
	std::string::size_type const s_length( s.length() );
	if ( w <= 0 ) {
		return std::string();
	} else if ( static_cast< int >( s_length ) > w ) {
		return s( 1, w );
	} else if ( static_cast< int >( s_length ) == w ) {
		return s;
	} else { // s_length < w
		return std::string( w - s_length, ' ' ) + s;
	}
}


/// @brief Fstring Format (Width of Fstring)
std::string
A( Fstring const & s )
{
	return s;
}


/// @brief Blank string
std::string
X( int const w )
{
	return std::string( std::max( w, 0 ), ' ' );
}


/// @brief Blank string
std::string
space( int const w )
{
	return std::string( std::max( w, 0 ), ' ' );
}

/// @brief Blank string
std::string
repeat( int const w, char c  )
{
	return std::string( std::max( w, 0 ), c );
}


// Logical Formatting


/// @brief Logical Format
std::string
L( int const w, bool const & t )
{
	if ( w <= 1 ) {
		return std::string( 1, ( t ? 'T' : 'F' ) );
	} else {
		return std::string( w-1, ' ' ) + std::string( 1, ( t ? 'T' : 'F' ) );
	}
}


/// @brief Logical Format (Width==1)
std::string
L( bool const & t )
{
	return std::string( 1, ( t ? 'T' : 'F' ) );
}


// Floating Point Formatting


/// @brief Exponential Format: float
std::string
E( int const w, int const d, float const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) ) << std::setw( w ) << t;
	return fmt_stream.str();
}


/// @brief Exponential Format: double
std::string
E( int const w, int const d, double const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) ) << std::setw( w ) << t;
	return fmt_stream.str();
}


/// @brief Exponential Format: long double
std::string
E( int const w, int const d, long double const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) ) << std::setw( w ) << t;
	return fmt_stream.str();
}


/// @brief Exponential Format: complex< float >
std::string
E( int const w, int const d, std::complex< float > const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) )
	 << '(' << std::setw( w ) << t.real() << ',' << std::setw( w ) << t.imag() << ')';
	return fmt_stream.str();
}


/// @brief Exponential Format: complex< double >
std::string
E( int const w, int const d, std::complex< double > const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) )
	 << '(' << std::setw( w ) << t.real() << ',' << std::setw( w ) << t.imag() << ')';
	return fmt_stream.str();
}


/// @brief Exponential Format: complex< long double >
std::string
E( int const w, int const d, std::complex< long double > const & t )
{
	if ( w <= 0 ) return std::string();
	std::ostringstream fmt_stream;
	fmt_stream << std::scientific << std::showpoint << std::uppercase
	 << std::setprecision( std::max( std::min( d, w-7 ), 0 ) )
	 << '(' << std::setw( w ) << t.real() << ',' << std::setw( w ) << t.imag() << ')';
	return fmt_stream.str();
}


/// @brief Fixed Point Format: float
std::string
F( int const w, int const d, float const & t )
{
	if ( w <= 0 ) return std::string();
	int const p( w - 3 + ( t >= 0.0f ? 1 : 0 ) + ( std::abs( t ) < 1.0f - 0.5f/std::pow( 10.0f, std::max( d, 0 ) ) ? 1 : 0 ) );
	std::stringstream fmt_stream;
	fmt_stream << std::fixed << std::showpoint
	 << std::setprecision( std::max( std::min( d, p ), 0 ) ) << std::setw( w ) << t;
	if ( ( t < 0.0f ) && ( t >= -0.5f ) ) { // Remove sign from -0.0
		float x;
		fmt_stream >> x;
		if ( x == 0.0f ) return F( w, d, 0.0f );
	}
	if ( fmt_stream.str().length() > std::string::size_type( w ) ) { // Exceeded field width
		if ( fmt_stream.str()[ 0 ] == '0' ) {
			return fmt_stream.str().substr( 1 ); // Trim lead zero
		} else if ( fmt_stream.str().substr( 0, 2 ) == "-0" ) {
			return '-' + fmt_stream.str().substr( 2 ); // Trim lead zero
		}
	}
	return fmt_stream.str();
}


/// @brief Fixed Point Format: double
std::string
F( int const w, int const d, double const & t )
{
	if ( w <= 0 ) return std::string();
	int const p( w - 3 + ( t >= 0.0 ? 1 : 0 ) + ( std::abs( t ) < 1.0 - 0.5/std::pow( 10.0, std::max( d, 0 ) ) ? 1 : 0 ) );
	std::stringstream fmt_stream;
	fmt_stream << std::fixed << std::showpoint
	 << std::setprecision( std::max( std::min( d, p ), 0 ) ) << std::setw( w ) << t;
	if ( ( t < 0.0 ) && ( t >= -0.5 ) ) { // Remove sign from -0.0
		double x;
		fmt_stream >> x;
		if ( x == 0.0 ) return F( w, d, 0.0 );
	}
	if ( fmt_stream.str().length() > std::string::size_type( w ) ) { // Exceeded field width
		if ( fmt_stream.str()[ 0 ] == '0' ) {
			return fmt_stream.str().substr( 1 ); // Trim lead zero
		} else if ( fmt_stream.str().substr( 0, 2 ) == "-0" ) {
			return '-' + fmt_stream.str().substr( 2 ); // Trim lead zero
		}
	}
	return fmt_stream.str();
}


/// @brief Fixed Point Format: long double
std::string
F( int const w, int const d, long double const & t )
{
	if ( w <= 0 ) return std::string();
	int const p( w - 3 + ( t >= 0.0l ? 1 : 0 ) + ( std::abs( t ) < 1.0l - 0.5l/std::pow( 10.0l, std::max( d, 0 ) ) ? 1 : 0 ) );
	std::stringstream fmt_stream;
	fmt_stream << std::fixed << std::showpoint
	 << std::setprecision( std::max( std::min( d, p ), 0 ) ) << std::setw( w ) << t;
	if ( ( t < 0.0L ) && ( t >= -0.5L ) ) { // Remove sign from -0.0
		long double x;
		fmt_stream >> x;
		if ( x == 0.0L ) return F( w, d, 0.0L );
	}
	if ( fmt_stream.str().length() > std::string::size_type( w ) ) { // Exceeded field width
		if ( fmt_stream.str()[ 0 ] == '0' ) {
			return fmt_stream.str().substr( 1 ); // Trim lead zero
		} else if ( fmt_stream.str().substr( 0, 2 ) == "-0" ) {
			return '-' + fmt_stream.str().substr( 2 ); // Trim lead zero
		}
	}
	return fmt_stream.str();
}


/// @brief Fixed Point Format: complex< float >
std::string
F( int const w, int const d, std::complex< float > const & t )
{
	return '(' + F( w, d, t.real() ) + ',' + F( w, d, t.imag() ) + ')';
}


/// @brief Fixed Point Format: complex< double >
std::string
F( int const w, int const d, std::complex< double > const & t )
{
	return '(' + F( w, d, t.real() ) + ',' + F( w, d, t.imag() ) + ')';
}


/// @brief Fixed Point Format: complex< long double >
std::string
F( int const w, int const d, std::complex< long double > const & t )
{
	return '(' + F( w, d, t.real() ) + ',' + F( w, d, t.imag() ) + ')';
}


/// @brief General Format: float
std::string
G( int const w, int const d, float const & t )
{
	double const m( std::abs( t ) );
	if ( m == 0.0 ) {
		return F( w-4, d-1, t ) + std::string( 4, ' ' );
	} else if ( m < 0.1 - ( 0.5 * std::pow( 10.0, -d-2 ) ) ) {
		return E( w, d, t );
	} else if ( m >= std::pow( 10.0, d ) - 0.5 ) {
		return E( w, d, t );
	} else {
		int const i( static_cast< int >( std::floor( std::log10( m / ( 0.1 - 0.5 / std::pow( 10.0, d+1 ) ) ) ) ) );
		return F( w-4, d-i, t ) + std::string( 4, ' ' );
	}
}


/// @brief General Format: double
std::string
G( int const w, int const d, double const & t )
{
	double const m( std::abs( t ) );
	if ( m == 0.0 ) {
		return F( w-4, d-1, t ) + std::string( 4, ' ' );
	} else if ( m < 0.1 - ( 0.5 * std::pow( 10.0, -d-2 ) ) ) {
		return E( w, d, t );
	} else if ( m >= std::pow( 10.0, d ) - 0.5 ) {
		return E( w, d, t );
	} else {
		int const i( static_cast< int >( std::floor( std::log10( m / ( 0.1 - 0.5 / std::pow( 10.0, d+1 ) ) ) ) ) );
		return F( w-4, d-i, t ) + std::string( 4, ' ' );
	}
}


/// @brief General Format: long double
std::string
G( int const w, int const d, long double const & t )
{
	long double const m( std::abs( t ) );
	if ( m == 0.0L ) {
		return F( w-4, d-1, t ) + std::string( 4, ' ' );
	} else if ( m < 0.1L - ( 0.5L * std::pow( 10.0L, -d-2 ) ) ) {
		return E( w, d, t );
	} else if ( m >= std::pow( 10.0L, d ) - 0.5L ) {
		return E( w, d, t );
	} else {
		int const i( static_cast< int >( std::floor( std::log10( m / ( 0.1l - 0.5l / std::pow( 10.0l, d+1 ) ) ) ) ) );
		return F( w-4, d-i, t ) + std::string( 4, ' ' );
	}
}


/// @brief General Format: complex< float >
std::string
G( int const w, int const d, std::complex< float > const & t )
{
	return '(' + G( w, d, t.real() ) + ',' + G( w, d, t.imag() ) + ')';
}


/// @brief General Format: complex< double >
std::string
G( int const w, int const d, std::complex< double > const & t )
{
	return '(' + G( w, d, t.real() ) + ',' + G( w, d, t.imag() ) + ')';
}


/// @brief General Format: complex< long double >
std::string
G( int const w, int const d, std::complex< long double > const & t )
{
	return '(' + G( w, d, t.real() ) + ',' + G( w, d, t.imag() ) + ')';
}


// Standard Formatting


/// @brief Standard Width Format: bool Specialization
template<>
std::string
SW( bool const & t )
{
	return L( 2, t );
}


/// @brief Standard Width Format: byte Specialization
template<>
std::string
SW( byte const & t )
{
	return I( 7, t );
}


/// @brief Standard Width Format: short Specialization
template<>
std::string
SW( short int const & t )
{
	return I( 7, t );
}


/// @brief Standard Width Format: unsigned short Specialization
template<>
std::string
SW( unsigned short int const & t )
{
	return I( 7, t );
}


/// @brief Standard Width Format: int Specialization
template<>
std::string
SW( int const & t )
{
	return I( 12, t );
}


/// @brief Standard Width Format: unsigned int Specialization
template<>
std::string
SW( unsigned int const & t )
{
	return I( 12, t );
}


/// @brief Standard Width Format: long int Specialization
template<>
std::string
SW( long int const & t )
{
	return I( 22, t );
}


/// @brief Standard Width Format: unsigned long int Specialization
template<>
std::string
SW( unsigned long int const & t )
{
	return I( 22, t );
}


/// @brief Standard Width Format: float Specialization
template<>
std::string
SW( float const & t )
{
	return G( 16, 7, t );
}


/// @brief Standard Width Format: double Specialization
template<>
std::string
SW( double const & t )
{
	return G( 26, 17, t );
}


/// @brief Standard Width Format: long double Specialization
template<>
std::string
SW( long double const & t )
{
	return G( 40, 32, t ); // This should be enough for 128-bit quad precision: long double is 80-bit on some platforms
}


/// @brief Standard Width Format: complex< float > Specialization
template<>
std::string
SW( std::complex< float > const & t )
{
	return '(' + SW( t.real() ) + ',' + SW( t.imag() ) + ')';
}


/// @brief Standard Width Format: complex< double > Specialization
template<>
std::string
SW( std::complex< double > const & t )
{
	return '(' + SW( t.real() ) + ',' + SW( t.imag() ) + ')';
}


/// @brief Standard Width Format: complex< long double > Specialization
template<>
std::string
SW( std::complex< long double > const & t )
{
	return '(' + SW( t.real() ) + ',' + SW( t.imag() ) + ')';
}


// Manipulators


/// @brief general: Manipulator to Turn Off scientific or fixed
std::ios_base &
general( std::ios_base & base )
{
	base.unsetf( std::ios_base::fixed );
	base.unsetf( std::ios_base::scientific );
	return base;
}


} // namespace format
} // namespace ObjexxFCL
