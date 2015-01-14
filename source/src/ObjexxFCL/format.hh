#ifndef INCLUDED_ObjexxFCL_format_hh
#define INCLUDED_ObjexxFCL_format_hh


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
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <sstream>
#include <string>


namespace ObjexxFCL {
namespace format {


// Constants
char const SPACE( ' ' );


// Types
typedef  char       *    cstring;
typedef  char const *  c_cstring;


// Formatted Input


// Bite: Inputs a String of Given Width from the Input Stream into a Value
template< typename T >
class Bite
{


public: // Creation


	/// @brief Width + Value Constructor
	inline
	Bite( int const w, T & t ) :
		w_( w ),
		d_( 0 ),
		t_( t )
	{}


	/// @brief Width + Precision + Value Constructor
	inline
	Bite( int const w, int const d, T & t ) :
		w_( w ),
		d_( d ),	// In Fortran d>=0 but that is not enforced here
		t_( t )
	{}


	/// @brief Destructor
	inline
	~Bite()
	{}


public: // I/O


	/// @brief Input a Bite from Stream
	friend
	inline
	std::istream &
	operator >>( std::istream & stream, Bite const & bite )
	{
		std::stringstream ss;
		char c;
		int i( 0 );
		while ( ( i < bite.w_ ) && ( stream ) && ( stream.peek() != '\n' ) ) {
			stream.get( c );
			if ( stream ) ss << c;
			++i;
		}
		bite.assign( ss );
		stream.setstate( stream.rdstate() | ( ss.rdstate() & ~std::ios_base::eofbit ) );
		return stream;
	}


private: // I/O


	/// @brief Assign Stream Bite to Value: Generic Implementation
	inline
	void
	assign( std::stringstream & ss ) const
	{
		ss >> t_;
	}


private: // Data


	int w_; // Width
	int d_; // Precision (controls implied decimal point location)
	T & t_; // Reference to value


}; // Bite


/// @brief Input a Bite from Stream
#ifndef __clang__
template< typename T >
std::istream &
operator >>( std::istream & stream, Bite< T > const & bite );
#endif


/// @brief string is Blank?
inline
bool
is_blank_string( std::string const & s )
{
	if ( s.empty() ) {
		return true;
	} else if ( s.find_first_not_of( ' ' ) == std::string::npos ) {
		return true;
	} else {
		return false;
	}
}


// Bite Explicit Specializations


	/// @brief Assign Stream Bite to Value: bool Specialization
	template<>
	inline
	void
	Bite< bool >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = false;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = false;
		}
	}


	/// @brief Assign Stream Bite to Value: byte Specialization
	template<>
	void
	Bite< byte >::assign( std::stringstream & ss ) const;


	/// @brief Assign Stream Bite to Value: ubyte Specialization
	template<>
	void
	Bite< ubyte >::assign( std::stringstream & ss ) const;


	/// @brief Assign Stream Bite to Value: char Specialization
	template<>
	inline
	void
	Bite< char >::assign( std::stringstream & ss ) const
	{
		t_ = ( is_blank_string( ss.str() ) ? SPACE : ss.str()[ 0 ] );
	}


	/// @brief Assign Stream Bite to Value: signed char Specialization
	template<>
	inline
	void
	Bite< signed char >::assign( std::stringstream & ss ) const
	{
		t_ = ( is_blank_string( ss.str() ) ? SPACE : ss.str()[ 0 ] );
	}


	/// @brief Assign Stream Bite to Value: unsigned char Specialization
	template<>
	inline
	void
	Bite< unsigned char >::assign( std::stringstream & ss ) const
	{
		t_ = ( is_blank_string( ss.str() ) ? SPACE : ss.str()[ 0 ] );
	}


	/// @brief Assign Stream Bite to Value: short int Specialization
	template<>
	inline
	void
	Bite< short int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: unsigned short int Specialization
	template<>
	inline
	void
	Bite< unsigned short int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: int Specialization
	template<>
	inline
	void
	Bite< int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: unsigned int Specialization
	template<>
	inline
	void
	Bite< unsigned int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: long int Specialization
	template<>
	inline
	void
	Bite< long int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: unsigned long int Specialization
	template<>
	inline
	void
	Bite< unsigned long int >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0;
		} else {
			ss >> t_;
			if ( ss.fail() ) t_ = 0;
		}
	}


	/// @brief Assign Stream Bite to Value: float Specialization
	template<>
	inline
	void
	Bite< float >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0f;
		} else {
			using std::string;
			string s( stripped( ss.str() ) );
			if ( s.find_first_of( "EeDd" ) == string::npos ) {	// No C++ supported exponent
				string::size_type const i( s.find_last_of( "+-" ) );
				if ( ( i != string::npos ) && ( i > 0 ) ) { // Fortran exponent form: num+ee or num-ee
					s.insert( i, 1, 'E' ); // Insert the E
					ss.str( s );
				}
			}
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0f;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0, static_cast< double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: double Specialization
	template<>
	inline
	void
	Bite< double >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0;
		} else {
			using std::string;
			string s( stripped( ss.str() ) );
			if ( s.find_first_of( "EeDd" ) == string::npos ) {	// No C++ supported exponent
				string::size_type const i( s.find_last_of( "+-" ) );
				if ( ( i != string::npos ) && ( i > 0 ) ) { // Fortran exponent form: num+ee or num-ee
					s.insert( i, 1, 'E' ); // Insert the E
					ss.str( s );
				}
			}
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0, static_cast< double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: long double Specialization
	template<>
	inline
	void
	Bite< long double >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0l;
		} else {
			using std::string;
			string s( stripped( ss.str() ) );
			if ( s.find_first_of( "EeDd" ) == string::npos ) {	// No C++ supported exponent
				string::size_type const i( s.find_last_of( "+-" ) );
				if ( ( i != string::npos ) && ( i > 0 ) ) { // Fortran exponent form: num+ee or num-ee
					s.insert( i, 1, 'E' ); // Insert the E
					ss.str( s );
				}
			}
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0l;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0l, static_cast< long double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: complex< float > Specialization
	template<>
	inline
	void
	Bite< std::complex< float > >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0f;
		} else {
			//Do Add support for the Fortran exponent forms: num+ee and num-ee
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0f;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == std::string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0, static_cast< double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: complex< double > Specialization
	template<>
	inline
	void
	Bite< std::complex< double > >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0;
		} else {
			//Do Add support for the Fortran exponent forms: num+ee and num-ee
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == std::string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0, static_cast< double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: complex< long double > Specialization
	template<>
	inline
	void
	Bite< std::complex< long double > >::assign( std::stringstream & ss ) const
	{
		if ( is_blank_string( ss.str() ) ) {
			t_ = 0.0l;
		} else {
			//Do Add support for the Fortran exponent forms: num+ee and num-ee
			ss >> t_;
			if ( ss.fail() ) {
				t_ = 0.0l;
			} else if ( ( d_ != 0 ) && ( ss.str().find( '.' ) == std::string::npos ) ) { // Implicit decimal point
				t_ *= std::pow( 10.0l, static_cast< long double >( -d_ ) );
			}
		}
	}


	/// @brief Assign Stream Bite to Value: string Specialization
	template<>
	inline
	void
	Bite< std::string >::assign( std::stringstream & ss ) const
	{
		t_ = ss.str();
	}


	/// @brief Assign Stream Bite to Value: Fstring Specialization
	template<>
	void
	Bite< Fstring >::assign( std::stringstream & ss ) const;


// Bite Makers


/// @brief Width + Value Bite Maker
template< typename T >
inline
Bite< T >
bite( int const w, T & t )
{
	return Bite< T >( w, t );
}


/// @brief Width + Precision + Value Bite Maker
template< typename T >
inline
Bite< T >
bite( int const w, int const d, T & t )
{
	return Bite< T >( w, d, t );
}


/// @brief bool Bite Maker: Take One Character
inline
Bite< bool >
bite( bool & t )
{
	return Bite< bool >( 1, t );
}


/// @brief char Bite Maker: Take One Character
inline
Bite< char >
bite( char & t )
{
	return Bite< char >( 1, t );
}


/// @brief string Bite Maker: Take Rest of Line
inline
Bite< std::string >
bite( std::string & t )
{
  return Bite< std::string >( (std::numeric_limits< int >::max)(), t );
}


/// @brief Fstring Bite Maker: Take Length of Fstring
Bite< Fstring >
bite( Fstring & t );


// Skip: Skips Over a Bite of Specified Width from the Input Stream
class Skip
{


public: // Creation


	/// @brief Constructor
	inline
	explicit
	Skip( int const w = 1 ) :
		w_( w )
	{}


	/// @brief Destructor
	inline
	~Skip()
	{}


public: // I/O


	/// @brief Input a Skip from Stream
	friend
	inline
	std::istream &
	operator >>( std::istream & stream, Skip const & skip )
	{
		char c;
		int i( 0 );
		while ( ( i < skip.w_ ) && ( stream ) && ( stream.peek() != '\n' ) ) {
			stream.get( c );
			++i;
		}
		return stream;
	}


private: // Data


	int w_; // Width


}; // Skip


/// @brief Input a Skip from Stream
#ifndef __clang__
std::istream &
operator >>( std::istream & stream, Skip const & skip );
#endif


// Skip Maker and Manipulator


/// @brief Skip Maker
inline
Skip
skip( int const w = 1 )
{
	return Skip( w );
}


/// @brief Skip Rest of Line and Line Terminator (Manipulator)
inline
std::istream &
skip( std::istream & stream )
{
  return stream.ignore( (std::numeric_limits< std::streamsize >::max)(), '\n' );
}


// Formatted Output


// General Formatting


/// @brief Single-Spaced Format
template< typename T >
inline
std::string
SS( T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << ' ' << t;
	return fmt_stream.str();
}


/// @brief Single-Spaced Format: bool Specialization
template<>
std::string
SS( bool const & t );


/// @brief Single-Spaced Format: float Specialization
template<>
std::string
SS( float const & t );


/// @brief Single-Spaced Format: double Specialization
template<>
std::string
SS( double const & t );


/// @brief Single-Spaced Format: long double Specialization
template<>
std::string
SS( long double const & t );


/// @brief Single-Spaced Format: complex< float > Specialization
template<>
std::string
SS( std::complex< float > const & t );


/// @brief Single-Spaced Format: complex< double > Specialization
template<>
std::string
SS( std::complex< double > const & t );


/// @brief Single-Spaced Format: complex< long double > Specialization
template<>
std::string
SS( std::complex< long double > const & t );


/// @brief Left-Justified Width-Specified Format
template< typename T >
inline
std::string
LJ( int const w, T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::setw( w ) << t;
	return fmt_stream.str();
}


/// @brief Right-Justified Width-Specified Format
template< typename T >
inline
std::string
RJ( int const w, T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::right << std::setw( w ) << t;
	return fmt_stream.str();
}


// String Formatting


/// @brief char Format
std::string
A( int const w, char const c );


/// @brief char Format (Width==1)
std::string
A( char const c );


/// @brief cstring Format
std::string
A( int const w, c_cstring const s );


/// @brief cstring Format (Width of cstring)
std::string
A( c_cstring const s );


/// @brief string Format
std::string
A( int const w, std::string const & s );


/// @brief string Format (Width of string)
std::string const &
A( std::string const & s );


/// @brief Fstring Format
std::string
A( int const w, Fstring const & s );


/// @brief Fstring Format (Width of Fstring)
std::string
A( Fstring const & s );


/// @brief Blank string
std::string
X( int const w );


/// @brief Blank string
std::string
space( int const w );

/// @brief w*c
std::string
repeat( int const w, char c );


// Logical Formatting


/// @brief Logical Format
std::string
L( int const w, bool const & t );


/// @brief Logical Format (Width==1)
std::string
L( bool const & t );


// Integer Formatting


/// @brief Integer Format
template< typename T >
inline
std::string
I( int const w, T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::right << std::setw( w ) << t;
	return fmt_stream.str();
}


/// @brief Integer Format with Minimum Digits
template< typename T >
inline
std::string
I( int const w, int const m, T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::right << std::setfill( '0' ) << std::setw( std::min( m, w ) ) << t;
	std::string const str( fmt_stream.str() );
	return std::string( std::max( w - static_cast< int >( str.length() ), 0 ), ' ' ) + str;
}


// Floating Point Formatting


/// @brief Exponential Format: float
std::string
E( int const w, int const d, float const & t );


/// @brief Exponential Format: double
std::string
E( int const w, int const d, double const & t );


/// @brief Exponential Format: long double
std::string
E( int const w, int const d, long double const & t );


/// @brief Exponential Format: complex< float >
std::string
E( int const w, int const d, std::complex< float > const & t );


/// @brief Exponential Format: complex< double >
std::string
E( int const w, int const d, std::complex< double > const & t );


/// @brief Exponential Format: complex< long double >
std::string
E( int const w, int const d, std::complex< long double > const & t );


/// @brief Fixed Point Format: float
std::string
F( int const w, int const d, float const & t );


/// @brief Fixed Point Format: double
std::string
F( int const w, int const d, double const & t );


/// @brief Fixed Point Format: long double
std::string
F( int const w, int const d, long double const & t );


/// @brief Fixed Point Format: complex< float >
std::string
F( int const w, int const d, std::complex< float > const & t );


/// @brief Fixed Point Format: complex< double >
std::string
F( int const w, int const d, std::complex< double > const & t );


/// @brief Fixed Point Format: complex< long double >
std::string
F( int const w, int const d, std::complex< long double > const & t );


/// @brief General Format: float
std::string
G( int const w, int const d, float const & t );


/// @brief General Format: double
std::string
G( int const w, int const d, double const & t );


/// @brief General Format: long double
std::string
G( int const w, int const d, long double const & t );


/// @brief General Format: complex< float >
std::string
G( int const w, int const d, std::complex< float > const & t );


/// @brief General Format: complex< double >
std::string
G( int const w, int const d, std::complex< double > const & t );


/// @brief General Format: complex< long double >
std::string
G( int const w, int const d, std::complex< long double > const & t );


// Standard Formatting


/// @brief Standard Width Format: Default Implementation
template< typename T >
inline
std::string
SW( T const & t )
{
	std::ostringstream fmt_stream;
	fmt_stream << std::left << std::noshowpoint << std::uppercase << t;
	return fmt_stream.str();
}


/// @brief Standard Width Format: bool Specialization
template<>
std::string
SW( bool const & t );


/// @brief Standard Width Format: byte Specialization
template<>
std::string
SW( byte const & t );


/// @brief Standard Width Format: short Specialization
template<>
std::string
SW( short int const & t );


/// @brief Standard Width Format: unsigned short Specialization
template<>
std::string
SW( unsigned short int const & t );


/// @brief Standard Width Format: int Specialization
template<>
std::string
SW( int const & t );


/// @brief Standard Width Format: unsigned int Specialization
template<>
std::string
SW( unsigned int const & t );


/// @brief Standard Width Format: long int Specialization
template<>
std::string
SW( long int const & t );


/// @brief Standard Width Format: unsigned long int Specialization
template<>
std::string
SW( unsigned long int const & t );


/// @brief Standard Width Format: float Specialization
template<>
std::string
SW( float const & t );


/// @brief Standard Width Format: double Specialization
template<>
std::string
SW( double const & t );


/// @brief Standard Width Format: long double Specialization
template<>
std::string
SW( long double const & t );


/// @brief Standard Width Format: complex< float > Specialization
template<>
std::string
SW( std::complex< float > const & t );


/// @brief Standard Width Format: complex< double > Specialization
template<>
std::string
SW( std::complex< double > const & t );


/// @brief Standard Width Format: complex< long double > Specialization
template<>
std::string
SW( std::complex< long double > const & t );


// Manipulators


/// @brief general: Manipulator to Turn Off scientific or fixed
std::ios_base &
general( std::ios_base & base );


// Utility


/// @brief Newline utility for formatted output implied DO loop emulation
inline
std::string
nl_if( int const i, int const n )
{
	if ( ( i > 1 ) && ( ( i - 1 ) % n == 0 ) ) {
		return "\n";
	} else {
		return "";
	}
}


} // namespace format
} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_format_HH
