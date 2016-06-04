#ifndef INCLUDED_ObjexxFCL_Cstring_hh
#define INCLUDED_ObjexxFCL_Cstring_hh


// Cstring: C String Wrapper
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
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iosfwd>
#include <string>

#ifdef CXX11
#include <type_traits>  // for swap
#else
#include <algorithm>
#endif

namespace ObjexxFCL {


// Types
typedef  char       *    cstring;
typedef  char const *  c_cstring;


/// @brief Cstring: C String Wrapper
///
/// @remarks
///  @li A memory-managed C string (char*) wrapper for convenience when using a C-style interfaces
///  @li Explicit conversion from std::string
///  @li Implicit conversion from/to char* for argument passing
///  @li Automatic memory management
///  @li Invariant: Null-terminated upon return from any constructor or function
class Cstring
{


public: // Types


	// STL Style
	typedef  std::size_t  size_type;

	// C++ Style
	typedef  std::size_t  Size;


public: // Creation


	/// @brief Default Constructor
	inline
	Cstring() :
		str_( new char[ 1 ] )
	{
		str_[ 0 ] ='\0';
	}


	/// @brief Copy Constructor
	inline
	Cstring( Cstring const & s ) :
		str_( new char[ std::strlen( s.str_ ) + 1 ] )
	{
		std::memcpy( str_, s.str_, std::strlen( s.str_ ) + 1 );
	}


	/// @brief C string Constructor: Implicit Conversion
	inline
	Cstring( c_cstring const s ) :
		str_( new char[ std::strlen( s ) + 1 ] )
	{
		std::memcpy( str_, s, std::strlen( s ) + 1 );
	}


	/// @brief std::string Constructor
	inline
	explicit
	Cstring( std::string const & s ) :
		str_( new char[ s.length() + 1 ] )
	{
		size_type const len( s.length() );
		s.copy( str_, len );
		str_[ len ] = '\0';
	}


	/// @brief Cstring + Length Constructor
	inline
	Cstring(
		Cstring const & s,
		size_type const len
	) :
		str_( new char[ len + 1 ] )
	{
		assert( len <= s.length() );
		std::memcpy( str_, s.str_, len );
		str_[ len ] = '\0';
	}


	/// @brief C string + Length Constructor
	inline
	Cstring(
		c_cstring const s,
		size_type const len
	) :
		str_( new char[ len + 1 ] )
	{
		assert( len <= std::strlen( s ) );
		std::memcpy( str_, s, len );
		str_[ len ] = '\0';
	}


	/// @brief std::string + Length Constructor
	inline
	Cstring(
		std::string const & s,
		size_type const len
	) :
		str_( new char[ len + 1 ] )
	{
		assert( len <= s.length() );
		s.copy( str_, len );
		str_[ len ] = '\0';
	}


	/// @brief char Constructor
	inline
	explicit
	Cstring( char const c ) :
		str_( new char[ 2 ] )
	{
		str_[ 0 ] = c;
		str_[ 1 ] = '\0';
	}


	/// @brief Length Constructor
	inline
	explicit
	Cstring( size_type const len ) :
		str_( new char[ len + 1 ] )
	{
		std::memset( str_, ' ', len );
		str_[ len ] = '\0';
	}


	/// @brief Length Constructor
	inline
	explicit
	Cstring( int const len ) :
		str_( new char[ len + 1 ] )
	{
		std::memset( str_, ' ', len );
		str_[ len ] = '\0';
	}


	/// @brief Destructor
	inline
	virtual
	~Cstring()
	{
		delete[] str_;
	}


public: // Conversion


	/// @brief C string Conversion: Invalid after str_ is reallocated
	inline
	operator c_cstring() const
	{
		return str_;
	}


	/// @brief C string Conversion: Invalid after str_ is reallocated
	inline
	operator cstring()
	{
		return str_;
	}


public: // Assignment


	/// @brief Copy Assignment
	inline
	Cstring &
	operator =( Cstring const & s )
	{
		if ( this != &s ) {
			size_type const len( s.length() + 1 );
			delete[] str_; str_ = new char[ len ];
			std::memcpy( str_, s.str_, len );
		}
		return *this;
	}


	/// @brief cstring Assignment
	inline
	Cstring &
	operator =( c_cstring const s )
	{
		size_type const len( std::strlen( s ) + 1 );
		delete[] str_; str_ = new char[ len ];
		std::memmove( str_, s, len );
		return *this;
	}


	/// @brief std::string Assignment
	inline
	Cstring &
	operator =( std::string const & s )
	{
		size_type const len( s.length() );
		delete[] str_; str_ = new char[ len + 1 ];
		s.copy( str_, len );
		str_[ len ] = '\0';
		return *this;
	}


	/// @brief char Assignment
	inline
	Cstring &
	operator =( char const c )
	{
		delete[] str_; str_ = new char[ 2 ];
		str_[ 0 ] = c;
		str_[ 1 ] = '\0';
		return *this;
	}


	/// @brief Cstring Append
	inline
	Cstring &
	operator +=( Cstring const & s )
	{
		Cstring( *this + s ).swap( *this );
		return *this;
	}


	/// @brief cstring Append
	inline
	Cstring &
	operator +=( c_cstring const s )
	{
		Cstring( *this + s ).swap( *this );
		return *this;
	}


	/// @brief std::string Append
	inline
	Cstring &
	operator +=( std::string const & s )
	{
		Cstring( *this + s ).swap( *this );
		return *this;
	}


	/// @brief char Append
	inline
	Cstring &
	operator +=( char const c )
	{
		Cstring( *this + c ).swap( *this );
		return *this;
	}


public: // Predicate


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return ( std::strlen( str_ ) == 0 );
	}


	/// @brief Blank?
	inline
	bool
	is_blank() const
	{
		return ( len_trim() == 0 );
	}


	/// @brief Not blank?
	inline
	bool
	not_blank() const
	{
		return ( len_trim() != 0 );
	}


	/// @brief Has Any Character of a Cstring?
	bool
	has_any_of( Cstring const & s ) const;


	/// @brief Has Any Character of a cstring?
	bool
	has_any_of( c_cstring const & s ) const;


	/// @brief Has Any Character of a std::string?
	bool
	has_any_of( std::string const & s ) const;


	/// @brief Has a Character?
	bool
	has_any_of( char const c ) const;


	/// @brief Has a Character?
	bool
	has( char const c ) const;


public: // Inspector


	/// @brief Length
	inline
	size_type
	length() const
	{
		return std::strlen( str_ );
	}


	/// @brief Length
	inline
	size_type
	len() const
	{
		return std::strlen( str_ );
	}


	/// @brief Size
	inline
	size_type
	size() const
	{
		return std::strlen( str_ );
	}


	/// @brief Length Space-Trimmed
	size_type
	len_trim() const;


	/// @brief Length Whitespace-Trimmed
	size_type
	len_trim_whitespace() const;


	/// @brief Find First Occurrence of a Character
	size_type
	find( char const c ) const;


	/// @brief Find Last Occurrence of a Character
	size_type
	find_last( char const c ) const;


public: // Modifier


	/// @brief Lowercase
	Cstring &
	lowercase();


	/// @brief Uppercase
	Cstring &
	uppercase();


	/// @brief Left Justify
	Cstring &
	left_justify();


	/// @brief Right Justify
	Cstring &
	right_justify();


	/// @brief Trim Trailing Space
	inline
	Cstring &
	trim()
	{
		str_[ len_trim() ] = '\0';
		return *this;
	}


	/// @brief Trim Trailing Whitespace
	inline
	Cstring &
	trim_whitespace()
	{
		str_[ len_trim_whitespace() ] = '\0';
		return *this;
	}


	/// @brief Center
	Cstring &
	center();


	/// @brief Compress Out Whitespace
	Cstring &
	compress();


	/// @brief swap( Cstring )
	inline
	void
	swap( Cstring & s )
	{
		std::swap( str_, s.str_ );
	}


	/// @brief swap( Cstring, Cstring )
	friend
	inline
	void
	swap( Cstring & s, Cstring & t )
	{
		std::swap( s.str_, t.str_ );
	}


public: // Subscript


	/// @brief Cstring[ i ] const
	inline
	char
	operator []( size_type const i ) const
	{
		assert( i < std::strlen( str_ ) );
		return str_[ i ];
	}


	/// @brief Cstring[ i ]
	inline
	char &
	operator []( size_type const i )
	{
		assert( i < std::strlen( str_ ) );
		return str_[ i ];
	}


	/// @brief Cstring[ i ] const
	/// @note  Overload prevents ambiguity with built-in operator[] with int arguments
	inline
	char
	operator []( int const i ) const
	{
		assert( i >= 0 );
		assert( static_cast< size_type >( i ) < std::strlen( str_ ) );
		return str_[ i ];
	}


	/// @brief Cstring[ i ]
	/// @note  Overload prevents ambiguity with built-in operator[] with int arguments
	inline
	char &
	operator []( int const i )
	{
		assert( i >= 0 );
		assert( static_cast< size_type >( i ) < std::strlen( str_ ) );
		return str_[ i ];
	}


public: // Concatenation


	/// @brief Cstring + Cstring
	friend
	inline
	Cstring
	operator +( Cstring const & s, Cstring const & t )
	{
		size_type const s_len( s.length() );
		size_type const t_len( t.length() );
		Cstring u( s_len + t_len );
		std::memcpy( u.str_, s.str_, s_len );
		std::memcpy( u.str_ + s_len, t.str_, t_len + 1 );
		return u;
	}


	/// @brief Cstring + cstring
	friend
	inline
	Cstring
	operator +( Cstring const & s, c_cstring const t )
	{
		size_type const s_len( s.length() );
		size_type const t_len( std::strlen( t ) );
		Cstring u( s_len + t_len );
		std::memcpy( u.str_, s.str_, s_len );
		std::memcpy( u.str_ + s_len, t, t_len + 1 );
		return u;
	}


	/// @brief cstring + Cstring
	friend
	inline
	Cstring
	operator +( c_cstring const s, Cstring const & t )
	{
		size_type const s_len( std::strlen( s ) );
		size_type const t_len( t.length() );
		Cstring u( s_len + t_len );
		std::memcpy( u.str_, s, s_len );
		std::memcpy( u.str_ + s_len, t.str_, t_len + 1 );
		return u;
	}


	/// @brief Cstring + std::string
	friend
	inline
	Cstring
	operator +( Cstring const & s, std::string const & t )
	{
		size_type const s_len( s.length() );
		size_type const t_len( t.length() );
		Cstring u( s_len + t_len );
		std::memcpy( u.str_, s.str_, s_len );
		t.copy( u.str_ + s_len, t_len );
		return u;
	}


	/// @brief Cstring + char
	friend
	inline
	Cstring
	operator +( Cstring const & s, char const c )
	{
		size_type const s_len( s.length() );
		Cstring u( s_len + 1 );
		std::memcpy( u.str_, s.str_, s_len );
		u.str_[ s_len ] = c;
		return u;
	}


	/// @brief char + Cstring
	friend
	inline
	Cstring
	operator +( char const c, Cstring const & t )
	{
		size_type const t_len( t.length() );
		Cstring u( 1 + t_len );
		u.str_[ 0 ] = c;
		std::memcpy( u.str_ + 1, t.str_, t_len + 1 );
		return u;
	}


public: // Generator


	/// @brief Lowercased Copy
	inline
	Cstring
	lowercased() const
	{
		return Cstring( *this ).lowercase();
	}


	/// @brief Uppercased Copy
	inline
	Cstring
	uppercased() const
	{
		return Cstring( *this ).uppercase();
	}


	/// @brief Left-Justified Copy
	inline
	Cstring
	left_justified() const
	{
		return Cstring( *this ).left_justify();
	}


	/// @brief Right-Justified Copy
	inline
	Cstring
	right_justified() const
	{
		return Cstring( *this ).right_justify();
	}


	/// @brief Space-Trimmed Copy
	inline
	Cstring
	trimmed() const
	{
		return Cstring( *this, len_trim() );
	}


	/// @brief Whitespace-Trimmed Copy
	inline
	Cstring
	trimmed_whitespace() const
	{
		return Cstring( *this, len_trim_whitespace() );
	}


	/// @brief Centered Copy
	inline
	Cstring
	centered() const
	{
		return Cstring( *this ).center();
	}


	/// @brief Compressed Copy
	inline
	Cstring
	compressed() const
	{
		return Cstring( *this ).compress();
	}


public: // Comparison


	/// @brief Cstring == Cstring
	friend
	inline
	bool
	operator ==( Cstring const & s, Cstring const & t )
	{
		return ( std::strcmp( s.str_, t.str_ ) == 0 );
	}


	/// @brief Cstring != Cstring
	friend
	inline
	bool
	operator !=( Cstring const & s, Cstring const & t )
	{
		return ( std::strcmp( s.str_, t.str_ ) != 0 );
	}


	/// @brief Cstring == cstring
	friend
	inline
	bool
	operator ==( Cstring const & s, c_cstring const t )
	{
		return ( std::strcmp( s.str_, t ) == 0 );
	}


	/// @brief cstring == Cstring
	friend
	inline
	bool
	operator ==( c_cstring const t, Cstring const & s )
	{
		return ( std::strcmp( s.str_, t ) == 0 );
	}


	/// @brief Cstring != cstring
	friend
	inline
	bool
	operator !=( Cstring const & s, c_cstring const t )
	{
		return ( std::strcmp( s.str_, t ) != 0 );
	}


	/// @brief cstring != Cstring
	friend
	inline
	bool
	operator !=( c_cstring const t, Cstring const & s )
	{
		return ( std::strcmp( s.str_, t ) != 0 );
	}


	/// @brief Cstring == std::string
	friend
	inline
	bool
	operator ==( Cstring const & s, std::string const & t )
	{
		return ( s.str_ == t );
	}


	/// @brief std::string == Cstring
	friend
	inline
	bool
	operator ==( std::string const & t, Cstring const & s )
	{
		return ( t == s.str_ );
	}


	/// @brief Cstring != std::string
	friend
	inline
	bool
	operator !=( Cstring const & s, std::string const & t )
	{
		return ( s.str_ != t );
	}


	/// @brief std::string != Cstring
	friend
	inline
	bool
	operator !=( std::string const & t, Cstring const & s )
	{
		return ( t != s.str_ );
	}


	/// @brief Cstring == char
	friend
	inline
	bool
	operator ==( Cstring const & s, char const c )
	{
		return ( ( s.length() == 1 ) && ( s.str_[ 0 ] == c ) );
	}


	/// @brief char == Cstring
	friend
	inline
	bool
	operator ==( char const c, Cstring const & s )
	{
		return ( ( s.length() == 1 ) && ( s.str_[ 0 ] == c ) );
	}


	/// @brief Cstring != char
	friend
	inline
	bool
	operator !=( Cstring const & s, char const c )
	{
		return ( ( s.length() != 1 ) || ( s.str_[ 0 ] != c ) );
	}


	/// @brief char != Cstring
	friend
	inline
	bool
	operator !=( char const c, Cstring const & s )
	{
		return ( ( s.length() != 1 ) || ( s.str_[ 0 ] != c ) );
	}


	/// @brief Cstring == Cstring Case-Insensitively?
	friend
	bool
	equali( Cstring const & s, Cstring const & t );


	/// @brief Cstring == cstring Case-Insensitively?
	friend
	bool
	equali( Cstring const & s, c_cstring const t );


	/// @brief cstring == Cstring Case-Insensitively?
	friend
	bool
	equali( c_cstring const s, Cstring const & t );


	/// @brief Cstring == std::string Case-Insensitively?
	friend
	bool
	equali( Cstring const & s, std::string const & t );


	/// @brief std::string == Cstring Case-Insensitively?
	friend
	bool
	equali( std::string const & s, Cstring const & t );


	/// @brief Cstring == char Case-Insensitively?
	friend
	bool
	equali( Cstring const & s, char const c );


	/// @brief char == Cstring Case-Insensitively?
	friend
	bool
	equali( char const c, Cstring const & s );


public: // I/O


	/// @brief Output to Stream
	friend
	std::ostream &
	operator <<( std::ostream & stream, Cstring const & s );


	/// @brief Input from Stream
	friend
	std::istream &
	operator >>( std::istream & stream, Cstring & s );


public: // Data


	static size_type const npos = static_cast< size_type >( -1 );


private: // Data


	/// @brief String
	char * str_;


}; // Cstring


/// @brief swap( Cstring, Cstring );
void
swap( Cstring & s, Cstring & t );


/// @brief Cstring + Cstring
Cstring
operator +( Cstring const & s, Cstring const & t );


/// @brief Cstring + cstring
Cstring
operator +( Cstring const & s, c_cstring const t );


/// @brief cstring + Cstring
Cstring
operator +( c_cstring const s, Cstring const & t );


/// @brief Cstring + std::string
Cstring
operator +( Cstring const & s, std::string const & t );


/// @brief Cstring + char
Cstring
operator +( Cstring const & s, char const c );


/// @brief char + Cstring
Cstring
operator +( char const c, Cstring const & t );


/// @brief Cstring == Cstring
bool
operator ==( Cstring const & s, Cstring const & t );


/// @brief Cstring != Cstring
bool
operator !=( Cstring const & s, Cstring const & t );


/// @brief Cstring == cstring
bool
operator ==( Cstring const & s, c_cstring const t );


/// @brief cstring == Cstring
bool
operator ==( c_cstring const t, Cstring const & s );


/// @brief Cstring != cstring
bool
operator !=( Cstring const & s, c_cstring const t );


/// @brief cstring != Cstring
bool
operator !=( c_cstring const t, Cstring const & s );


/// @brief Cstring == std::string
bool
operator ==( Cstring const & s, std::string const & t );


/// @brief std::string == Cstring
bool
operator ==( std::string const & t, Cstring const & s );


/// @brief Cstring != std::string
bool
operator !=( Cstring const & s, std::string const & t );


/// @brief std::string != Cstring
bool
operator !=( std::string const & t, Cstring const & s );


/// @brief Cstring == char
bool
operator ==( Cstring const & s, char const c );


/// @brief char == Cstring
bool
operator ==( char const c, Cstring const & s );


/// @brief Cstring != char
bool
operator !=( Cstring const & s, char const c );


/// @brief char != Cstring
bool
operator !=( char const c, Cstring const & s );


/// @brief Cstring == Cstring Case-Insensitively?
bool
equali( Cstring const & s, Cstring const & t );


/// @brief Cstring == cstring Case-Insensitively?
bool
equali( Cstring const & s, c_cstring const t );


/// @brief cstring == Cstring Case-Insensitively?
bool
equali( c_cstring const s, Cstring const & t );


/// @brief Cstring == std::string Case-Insensitively?
bool
equali( Cstring const & s, std::string const & t );


/// @brief std::string == Cstring Case-Insensitively?
bool
equali( std::string const & s, Cstring const & t );


/// @brief Cstring == char Case-Insensitively?
bool
equali( Cstring const & s, char const c );


/// @brief char == Cstring Case-Insensitively?
bool
equali( char const c, Cstring const & s );


/// @brief Output to Stream
std::ostream &
operator <<( std::ostream & stream, Cstring const & s );


/// @brief Input from Stream
std::istream &
operator >>( std::istream & stream, Cstring & s );


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_Cstring_HH
