#ifndef INCLUDED_ObjexxFCL_Fstring_hh
#define INCLUDED_ObjexxFCL_Fstring_hh


// Fstring: Fixed-Length Fortran-Compatible String and Substring
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
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/TypeTraits.hh>

// C++ Headers
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <iomanip>
#include <iosfwd>
#include <sstream>
#include <string>


namespace ObjexxFCL {


// Types
typedef  char       *    cstring;
typedef  char const *  c_cstring;


/// @brief Fstring: Fixed-Length Fortran-Compatible String
///
/// @remarks
///  @li Subscripts run from 1 to the length
///  @li Space-padding is used in comparisons and assignments
///  @li Internal string rep is not null-terminated
///  @li Zero-length Fstrings are supported but cannot be indexed into (no valid indices)
///  @li All the length constructors are needed to avoid ambiguity with the char constructor
///  @li Assignment can set length/string if Fstring is uninitialized (default constructed)
///  @li Substrings: Use s( i, j ) or s( i ) / Pass s( i, j ).ref() to a non-const Fstring& argument
///  @li Assumes that char is a single-byte ASCII-collated character
class Fstring
{


private: // Friend


	friend class Fsubstring;


public: // Types

	// STL Style
	typedef  std::size_t  size_type;
	typedef  void (*initializer_function)( Fstring & );

	// C++ Style
	typedef  std::size_t  Size;
	typedef  void (*InitializerFunction)( Fstring & );


public: // Creation


	/// @brief Default Constructor
	inline
	Fstring() :
		len_( 0 ),
		str_( 0 ),
		c_str_( 0 ),
		sub_( false )
	{}


	/// @brief Copy Constructor
	Fstring( Fstring const & s );


	/// @brief string Constructor
	Fstring( std::string const & s );


	/// @brief cstring Constructor
	Fstring( c_cstring const s );


	/// @brief char Constructor
	inline
	explicit
	Fstring( char const c ) :
		len_( 1 ),
		str_( new char[ 1 ] ),
		c_str_( 0 ),
		sub_( false )
	{
		str_[ 0 ] = c;
	}


	/// @brief signed char Constructor
	inline
	explicit
	Fstring( signed char const c ) :
		len_( 1 ),
		str_( new char[ 1 ] ),
		c_str_( 0 ),
		sub_( false )
	{
		str_[ 0 ] = static_cast< char >( c );
	}


	/// @brief unsigned char Constructor
	inline
	explicit
	Fstring( unsigned char const c ) :
		len_( 1 ),
		str_( new char[ 1 ] ),
		c_str_( 0 ),
		sub_( false )
	{
		str_[ 0 ] = static_cast< char >( c );
	}


	/// @brief Length Constructor
	explicit
	Fstring( short int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( long int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( unsigned short int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( unsigned int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( unsigned long int const len_a );


	/// @brief Length Constructor
	explicit
	Fstring( unsigned long long const len_a );


	/// @brief Length + Fstring Constructor
	Fstring( size_type const len_a, Fstring const & s );


	/// @brief Length + string Constructor
	Fstring( size_type const len_a, std::string const & s );


	/// @brief Length + cstring Constructor
	Fstring( size_type const len_a, c_cstring const s );


	/// @brief Length + char Constructor
	/// @note Fills with specified char => Use Fstring( len_a, "c" ) for space-padded single character
	Fstring( size_type const len_a, char const c );


	/// @brief Length + Initializer Constructor
	Fstring( size_type const len_a, initializer_function init );


	/// @brief Destructor
	inline
	virtual
	~Fstring()
	{
		if ( ! sub_ ) delete[] str_; // Substrings don't own/delete data
		delete[] c_str_;
	}


protected: // Creation


	/// @brief Substring Range Constructor
	Fstring( Fstring const & s, size_type const i, size_type const j );


	/// @brief Substring Tail Constructor
	Fstring( Fstring const & s, size_type const i );


public: // Conversion


	/// @brief string Conversion
	inline
	operator std::string() const
	{
		return std::string( str_, len_ );
	}


public: // Assignment


	/// @brief Copy Assignment
	Fstring &
	operator =( Fstring const & s );


	/// @brief = string
	Fstring &
	operator =( std::string const & s );


	/// @brief = cstring
	Fstring &
	operator =( c_cstring const s );


	/// @brief = char
	Fstring &
	operator =( char const c );


public: // Subscript


	/// @brief Constant char: s[ i ]
	inline
	char
	operator []( size_type const i ) const
	{
		assert( i > 0 );
		assert( i <= len_ );
		return str_[ i - 1 ];
	}


	/// @brief char: s[ i ]
	inline
	char &
	operator []( size_type const i )
	{
		assert( i > 0 );
		assert( i <= len_ );
		return str_[ i - 1 ];
	}


public: // Predicate


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return ( len_ == 0 );
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
		return ( len_trim() > 0 );
	}


	/// @brief Whitespace?
	inline
	bool
	is_whitespace() const
	{
		return ( len_trim_whitespace() == 0 );
	}


	/// @brief Not whitespace?
	inline
	bool
	not_whitespace() const
	{
		return ( len_trim_whitespace() > 0 );
	}


	/// @brief Has an Fstring?
	bool
	has( Fstring const & s ) const;


	/// @brief Has a string?
	bool
	has( std::string const & s ) const;


	/// @brief Has a cstring?
	bool
	has( c_cstring const s ) const;


	/// @brief Has a Character?
	bool
	has( char const c ) const;


	/// @brief Has Any Character of an Fstring?
	bool
	has_any_of( Fstring const & s ) const;


	/// @brief Has Any Character of a string?
	bool
	has_any_of( std::string const & s ) const;


	/// @brief Has Any Character of a cstring?
	bool
	has_any_of( c_cstring const s ) const;


	/// @brief Has a Character?
	bool
	has_any_of( char const c ) const;


	/// @brief Has a Prefix Case-Optionally?
	bool
	has_prefix( Fstring const & s, bool const exact_case = true ) const;


	/// @brief Has a Prefix Case-Optionally?
	bool
	has_prefix( c_cstring const s, bool const exact_case = true ) const;


	/// @brief Fstring is Readable as a Type Supporting Stream Input?
	template< typename T >
	inline
	bool
	is_type() const
	{
		if ( is_whitespace() ) { // Don't accept empty or whitespace Fstring
			return false;
		} else { // Try to read the Fstring as a T
			size_type b, e;
			size_type const trimmed_whitespace_length( trimmed_whitespace_range( b, e ) );
			std::istringstream t_stream( std::string( str_ + b - 1, trimmed_whitespace_length ) );
			T t;
			t_stream >> t;
			return ( ( t_stream ) && ( t_stream.eof() ) );
		}
	}


	/// @brief Fstring is Readable as a bool?
	inline
	bool
	is_bool() const
	{
		return is_type< bool >();
	}


	/// @brief Fstring is Readable as a short int?
	inline
	bool
	is_short() const
	{
		return is_type< short int >();
	}


	/// @brief Fstring is Readable as an int?
	inline
	bool
	is_int() const
	{
		return is_type< int >();
	}


	/// @brief Fstring is Readable as a long int?
	inline
	bool
	is_long() const
	{
		return is_type< long int >();
	}


	/// @brief Fstring is Readable as a unsigned short int?
	inline
	bool
	is_ushort() const
	{
		return is_type< unsigned short int >();
	}


	/// @brief Fstring is Readable as an unsigned int?
	inline
	bool
	is_uint() const
	{
		return is_type< unsigned int >();
	}


	/// @brief Fstring is Readable as a unsigned long int?
	inline
	bool
	is_ulong() const
	{
		return is_type< unsigned long int >();
	}


	/// @brief Fstring is Readable as a float?
	inline
	bool
	is_float() const
	{
		return is_type< float >();
	}


	/// @brief Fstring is Readable as a double?
	inline
	bool
	is_double() const
	{
		return is_type< double >();
	}


	/// @brief Fstring is Readable as a long double?
	inline
	bool
	is_longdouble() const
	{
		return is_type< long double >();
	}


	/// @brief Fstring is Readable as a char?
	inline
	bool
	is_char() const
	{
		return ( size() == 1 );
	}


	/// @brief Fstring is Readable as a string?
	inline
	bool
	is_string() const
	{
		return true;
	}


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return len_;
	}


	/// @brief Length
	inline
	size_type
	length() const
	{
		return len_;
	}


	/// @brief Length
	inline
	size_type
	len() const
	{
		return len_;
	}


	/// @brief Length Space-Trimmed
	size_type
	len_trim() const;


	/// @brief Length Whitespace-Trimmed
	size_type
	len_trim_whitespace() const;


	/// @brief Find First Occurrence of a Whitespace Character
	size_type
	find_whitespace() const;


	/// @brief Find First Occurrence of a Non-Whitespace Character
	size_type
	find_non_whitespace() const;


	/// @brief Find Last Occurrence of a Whitespace Character
	size_type
	find_last_whitespace() const;


	/// @brief Find Last Occurrence of a Non-Whitespace Character
	size_type
	find_last_non_whitespace() const;


	/// @brief Get Range of Whitespace-Trimmed Portion and Return its Length
	size_type
	trimmed_whitespace_range( size_type & b, size_type & e ) const;


	/// @brief Find First Occurrence of an Fstring
	size_type
	find( Fstring const & s ) const;


	/// @brief Find First Occurrence of a string
	size_type
	find( std::string const & s ) const;


	/// @brief Find First Occurrence of a cstring
	size_type
	find( c_cstring const s ) const;


	/// @brief Find First Occurrence of a Character
	size_type
	find( char const c ) const;


	/// @brief Find Last Occurrence of an Fstring
	size_type
	find_last( Fstring const & s ) const;


	/// @brief Find Last Occurrence of a string
	size_type
	find_last( std::string const & s ) const;


	/// @brief Find Last Occurrence of a cstring
	size_type
	find_last( c_cstring const s ) const;


	/// @brief Find Last Occurrence of a Character
	size_type
	find_last( char const c ) const;


	/// @brief Find First Occurrence of Any Character of an Fstring
	size_type
	find_first_of( Fstring const & s ) const;


	/// @brief Find First Occurrence of Any Character of a string
	size_type
	find_first_of( std::string const & s ) const;


	/// @brief Find First Occurrence of Any Character of a cstring
	size_type
	find_first_of( c_cstring const s ) const;


	/// @brief Find First Occurrence of a Character
	size_type
	find_first_of( char const c ) const;


	/// @brief Find First Occurrence of Any Character not of an Fstring
	size_type
	find_first_not_of( Fstring const & s ) const;


	/// @brief Find First Occurrence of Any Character not of a string
	size_type
	find_first_not_of( std::string const & s ) const;


	/// @brief Find First Occurrence of Any Character not of a cstring
	size_type
	find_first_not_of( c_cstring const s ) const;


	/// @brief Find First Occurrence of not a Character
	size_type
	find_first_not_of( char const c ) const;


	/// @brief Find Last Occurrence of Any Character of an Fstring
	size_type
	find_last_of( Fstring const & s ) const;


	/// @brief Find Last Occurrence of Any Character of a string
	size_type
	find_last_of( std::string const & s ) const;


	/// @brief Find Last Occurrence of Any Character of a cstring
	size_type
	find_last_of( c_cstring const s ) const;


	/// @brief Find Last Occurrence of a Character
	size_type
	find_last_of( char const c ) const;


	/// @brief Find Last Occurrence of Any Character not of an Fstring
	size_type
	find_last_not_of( Fstring const & s ) const;


	/// @brief Find Last Occurrence of Any Character not of a string
	size_type
	find_last_not_of( std::string const & s ) const;


	/// @brief Find Last Occurrence of Any Character not of a cstring
	size_type
	find_last_not_of( c_cstring const s ) const;


	/// @brief Find Last Occurrence not of a Character
	size_type
	find_last_not_of( char const c ) const;


	/// @brief Type of an Fstring for Type Supporting Stream Input
	template< typename T >
	inline
	T
	type_of() const
	{
		size_type b, e;
		size_type const trimmed_whitespace_length( trimmed_whitespace_range( b, e ) );
		std::istringstream t_stream( std::string( str_ + b - 1, trimmed_whitespace_length ) );
		T t;
		t_stream >> t;
		return ( ( t_stream ) && ( t_stream.eof() ) ? t : T() ); // Check is_type first
	}


	/// @brief short int of the Fstring
	inline
	short int
	short_of() const
	{
		return type_of< short int >();
	}


	/// @brief int of the Fstring
	inline
	int
	int_of() const
	{
		return type_of< int >();
	}


	/// @brief long int of the Fstring
	inline
	long int
	long_of() const
	{
		return type_of< long int >();
	}


	/// @brief unsigned short int of the Fstring
	inline
	unsigned short int
	ushort_of() const
	{
		return type_of< unsigned short int >();
	}


	/// @brief unsigned int of the Fstring
	inline
	unsigned int
	uint_of() const
	{
		return type_of< unsigned int >();
	}


	/// @brief unsigned long int of the Fstring
	inline
	unsigned long int
	ulong_of() const
	{
		return type_of< unsigned long int >();
	}


	/// @brief float of the Fstring
	inline
	float
	float_of() const
	{
		return type_of< float >();
	}


	/// @brief double of the Fstring
	inline
	double
	double_of() const
	{
		return type_of< double >();
	}


	/// @brief long double of the Fstring
	inline
	long double
	longdouble_of() const
	{
		return type_of< long double >();
	}


	/// @brief char of the Fstring
	inline
	char
	char_of() const
	{
		return ( len_ == 1 ? str_[ 0 ] : char() ); // Check is_type first
	}


	/// @brief string of the Fstring
	inline
	std::string
	string_of() const
	{
		return std::string( str_, len_ );
	}


public: // Modifier


	/// @brief Lowercase
	Fstring &
	lowercase();


	/// @brief Uppercase
	Fstring &
	uppercase();


	/// @brief Left Justify
	Fstring &
	left_justify();


	/// @brief Right Justify
	Fstring &
	right_justify();


	/// @brief Center
	Fstring &
	center();


	/// @brief Compress Out Whitespace
	Fstring &
	compress();


	/// @brief Trim Trailing Space
	/// @note  No effect for Fstring: Included for interface consistency
	inline
	Fstring &
	trim()
	{
		return *this;
	}


	/// @brief Trim Trailing Whitespace Replacing it with Space
	Fstring &
	trim_whitespace();


	/// @brief Strip Specified Characters from an Fstring's Tails
	Fstring &
	strip( std::string const & chars );


	/// @brief Strip Specified Characters from an Fstring's Left Tail
	Fstring &
	lstrip( std::string const & chars );


	/// @brief Strip Specified Characters from an Fstring's Right Tail
	Fstring &
	rstrip( std::string const & chars );


	/// @brief Strip Space from an Fstring's Tails
	Fstring &
	strip();


	/// @brief Strip Space from an Fstring's Left Tail
	Fstring &
	lstrip();


	/// @brief Strip Space from an Fstring's Right Tail
	Fstring &
	rstrip();


	/// @brief Strip Whitespace from an Fstring's Tails
	Fstring &
	strip_whitespace();


	/// @brief Strip Whitespace from an Fstring's Left Tail
	Fstring &
	lstrip_whitespace();


	/// @brief Strip Whitespace from an Fstring's Right Tail
	Fstring &
	rstrip_whitespace();


	/// @brief Clear
	inline
	Fstring &
	clear()
	{
		std::memset( str_, ' ', len_ );
		return *this;
	}


	/// @brief Overlay an Fstring
	Fstring &
	overlay( Fstring const & s, size_type const pos = 1 );


	/// @brief Overlay a string
	Fstring &
	overlay( std::string const & s, size_type const pos = 1 );


	/// @brief Overlay a cstring
	Fstring &
	overlay( c_cstring const s, size_type const pos = 1 );


public: // Generator


	/// @brief Left-Justified Copy
	inline
	Fstring
	left_justified() const
	{
		return Fstring( *this ).left_justify();
	}


	/// @brief Right-Justified Copy
	inline
	Fstring
	right_justified() const
	{
		return Fstring( *this ).right_justify();
	}


	/// @brief Centered Copy
	inline
	Fstring
	centered() const
	{
		return Fstring( *this ).center();
	}


	/// @brief Compressed Copy
	inline
	Fstring
	compressed() const
	{
		return Fstring( *this ).compress();
	}


	/// @brief Lowercased Copy
	inline
	Fstring
	lowercased() const
	{
		return Fstring( *this ).lowercase();
	}


	/// @brief Uppercased Copy
	inline
	Fstring
	uppercased() const
	{
		return Fstring( *this ).uppercase();
	}


	/// @brief Trailing Space Trimmed Copy
	inline
	Fstring
	trimmed() const
	{
		return Fstring( *this, 1, len_trim() );
	}


	/// @brief Trailing Whitespace Trimmed Copy
	inline
	Fstring
	trimmed_whitespace() const
	{
		return Fstring( *this, 1, len_trim_whitespace() );
	}


	/// @brief Specified Characters Stripped from Tails Copy
	inline
	Fstring
	stripped( std::string const & chars ) const
	{
		size_type const ib( find_first_not_of( chars ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib, find_last_not_of( chars ) );
		} else {
			return Fstring();
		}
	}


	/// @brief Specified Characters Stripped from Left Tail Copy
	inline
	Fstring
	lstripped( std::string const & chars ) const
	{
		size_type const ib( find_first_not_of( chars ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib );
		} else {
			return Fstring();
		}
	}


	/// @brief Specified Characters Stripped from Right Tail Copy
	inline
	Fstring
	rstripped( std::string const & chars ) const
	{
		return Fstring( *this, 1, find_last_not_of( chars ) );
	}


	/// @brief Space Stripped from Tails Copy
	inline
	Fstring
	stripped() const
	{
		size_type const ib( find_first_not_of( ' ' ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib, find_last_not_of( ' ' ) );
		} else {
			return Fstring();
		}
	}


	/// @brief Space Stripped from Left Tail Copy
	inline
	Fstring
	lstripped() const
	{
		size_type const ib( find_first_not_of( ' ' ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib );
		} else {
			return Fstring();
		}
	}


	/// @brief Space Stripped from Right Tail Copy
	inline
	Fstring
	rstripped() const
	{
		return Fstring( *this, 1, find_last_not_of( ' ' ) );
	}


	/// @brief Whitespace Stripped from Tails Copy
	inline
	Fstring
	stripped_whitespace() const
	{
		size_type const ib( find_first_not_of( " \t\000" ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib, find_last_not_of( " \t\000" ) );
		} else {
			return Fstring();
		}
	}


	/// @brief Whitespace Stripped from Left Tail Copy
	inline
	Fstring
	lstripped_whitespace() const
	{
		size_type const ib( find_first_not_of( " \t\000" ) );
		if ( ib > 0 ) {
			return Fstring( *this, ib );
		} else {
			return Fstring();
		}
	}


	/// @brief Whitespace Stripped from Right Tail Copy
	inline
	Fstring
	rstripped_whitespace() const
	{
		return Fstring( *this, 1, find_last_not_of( " \t\000" ) );
	}


	/// @brief Null-Terminated cstring Copy of the Fstring that is Owned by the Fstring
	c_cstring
	c_str() const;


	/// @brief Whitespace-Trimmed Null-Terminated cstring Copy of the Fstring that is Owned by the Fstring
	c_cstring
	t_str() const;


	/// @brief Non-Null-Terminated cstring Copy of the Fstring Data
	inline
	c_cstring
	data() const
	{
		return str_;
	}


	/// @brief Copy to a Pre-Allocated String
	size_type
	copy( cstring str, size_type const len_a, size_type const off = 0 ) const;


public: // Concatenation


	/// @brief Fstring + Fstring
	friend
	inline
	Fstring
	operator +( Fstring const & s, Fstring const & t )
	{
		Fstring u( static_cast< size_type >( s.len_ + t.len_ ) );
		std::memcpy( u.str_, s.str_, s.len_ );
		std::memcpy( u.str_ + s.len_, t.str_, t.len_ );
		return u;
	}


	/// @brief Fstring + string
	friend
	inline
	std::string
	operator +( Fstring const & s, std::string const & t )
	{
		return ( static_cast< std::string >( s ) + t );
	}


	/// @brief string + Fstring
	friend
	inline
	std::string
	operator +( std::string const & t, Fstring const & s )
	{
		return ( t + static_cast< std::string >( s ) );
	}


	/// @brief Fstring + cstring
	friend
	inline
	Fstring
	operator +( Fstring const & s, c_cstring const t )
	{
		size_type const t_len( std::strlen( t ) );
		Fstring u( s.len_ + t_len );
		std::memcpy( u.str_, s.str_, s.len_ );
		std::memcpy( u.str_ + s.len_, t, t_len );
		return u;
	}


	/// @brief cstring + Fstring
	friend
	inline
	Fstring
	operator +( c_cstring const s, Fstring const & t )
	{
		size_type const s_len( std::strlen( s ) );
		Fstring u( s_len + t.len_ );
		std::memcpy( u.str_, s, s_len );
		std::memcpy( u.str_ + s_len, t.str_, t.len_ );
		return u;
	}


	/// @brief Fstring + char
	friend
	inline
	Fstring
	operator +( Fstring const & s, char const c )
	{
		Fstring u( s.len_ + 1 );
		std::memcpy( u.str_, s.str_, s.len_ );
		u.str_[ s.len_ ] = c;
		return u;
	}


	/// @brief char + Fstring
	friend
	inline
	Fstring
	operator +( char const c, Fstring const & s )
	{
		Fstring u( 1 + s.len_ );
		u.str_[ 0 ] = c;
		std::memcpy( u.str_ + 1, s.str_, s.len_ );
		return u;
	}


public: // Comparison


	/// @brief Fstring == Fstring
	friend
	bool
	operator ==( Fstring const & s, Fstring const & t );


	/// @brief Fstring != Fstring
	friend
	inline
	bool
	operator !=( Fstring const & s, Fstring const & t )
	{
		return !( s == t );
	}


	/// @brief Fstring == string
	friend
	bool
	operator ==( Fstring const & s, std::string const & t );


	/// @brief string == Fstring
	friend
	inline
	bool
	operator ==( std::string const & t, Fstring const & s )
	{
		return ( s == t );
	}


	/// @brief Fstring != string
	friend
	inline
	bool
	operator !=( Fstring const & s, std::string const & t )
	{
		return !( s == t );
	}


	/// @brief string != Fstring
	friend
	inline
	bool
	operator !=( std::string const & t, Fstring const & s )
	{
		return !( s == t );
	}


	/// @brief Fstring == cstring
	friend
	bool
	operator ==( Fstring const & s, c_cstring const t );


	/// @brief cstring == Fstring
	friend
	inline
	bool
	operator ==( c_cstring const t, Fstring const & s )
	{
		return ( s == t );
	}


	/// @brief Fstring != cstring
	friend
	inline
	bool
	operator !=( Fstring const & s, c_cstring const t )
	{
		return !( s == t );
	}


	/// @brief cstring != Fstring
	friend
	inline
	bool
	operator !=( c_cstring const t, Fstring const & s )
	{
		return !( s == t );
	}


	/// @brief Fstring == char
	friend
	bool
	operator ==( Fstring const & s, char const c );


	/// @brief char == Fstring
	friend
	inline
	bool
	operator ==( char const c, Fstring const & s )
	{
		return ( s == c );
	}


	/// @brief Fstring != char
	friend
	inline
	bool
	operator !=( Fstring const & s, char const c )
	{
		return !( s == c );
	}


	/// @brief char != Fstring
	friend
	inline
	bool
	operator !=( char const c, Fstring const & s )
	{
		return !( s == c );
	}


	/// @brief Fstring == Fstring Case-Insensitively?
	friend
	inline
	bool
	equali( Fstring const & s, Fstring const & t )
	{
		return ( s.lowercased() == t.lowercased() );
	}


	/// @brief Fstring == string Case-Insensitively?
	friend
	inline
	bool
	equali( Fstring const & s, std::string const & t )
	{
		return ( s.lowercased() == ObjexxFCL::lowercased( t ) );
	}


	/// @brief string == Fstring Case-Insensitively?
	friend
	inline
	bool
	equali( std::string const & s, Fstring const & t )
	{
		return ( ObjexxFCL::lowercased( s ) == t.lowercased() );
	}


	/// @brief Fstring == char Case-Insensitively?
	friend
	inline
	bool
	equali( Fstring const & s, char const c )
	{
		return ( s.lowercased() == std::tolower( c ) );
	}


	/// @brief char == Fstring Case-Insensitively?
	friend
	inline
	bool
	equali( char const c, Fstring const & s )
	{
		return ( s.lowercased() == std::tolower( c ) );
	}


	/// @brief Fstring == Fstring Case-Optionally?
	inline
	bool
	equal( Fstring const & s, Fstring const & t, bool const exact_case = true )
	{
		if ( exact_case ) {
			return ( s == t );
		} else {
			return ( s.lowercased() == t.lowercased() );
		}
	}


	/// @brief Fstring == char Case-Optionally?
	inline
	bool
	equal( Fstring const & s, char const c, bool const exact_case = true )
	{
		if ( exact_case ) {
			return ( s == c );
		} else {
			return ( s.lowercased() == std::tolower( c ) );
		}
	}


	/// @brief char == Fstring Case-Optionally?
	inline
	bool
	equal( char const c, Fstring const & s, bool const exact_case = true )
	{
		if ( exact_case ) {
			return ( s == c );
		} else {
			return ( s.lowercased() == std::tolower( c ) );
		}
	}


	/// @brief Fstring <= Fstring
	friend
	bool
	operator <=( Fstring const & s, Fstring const & t );


	/// @brief Fstring < Fstring
	friend
	bool
	operator <( Fstring const & s, Fstring const & t );


	/// @brief Fstring >= Fstring
	friend
	inline
	bool
	operator >=( Fstring const & s, Fstring const & t )
	{
		return !( s < t );
	}


	/// @brief Fstring > Fstring
	friend
	inline
	bool
	operator >( Fstring const & s, Fstring const & t )
	{
		return !( s <= t );
	}


	/// @brief Fstring <= string
	friend
	bool
	operator <=( Fstring const & s, std::string const & t );


	/// @brief Fstring < string
	friend
	bool
	operator <( Fstring const & s, std::string const & t );


	/// @brief Fstring >= string
	friend
	inline
	bool
	operator >=( Fstring const & s, std::string const & t )
	{
		return !( s < t );
	}


	/// @brief Fstring > string
	friend
	inline
	bool
	operator >( Fstring const & s, std::string const & t )
	{
		return !( s <= t );
	}


	/// @brief string <= Fstring
	friend
	inline
	bool
	operator <=( std::string const & s, Fstring const & t )
	{
		return ( t >= s );
	}


	/// @brief string < Fstring
	friend
	inline
	bool
	operator <( std::string const & s, Fstring const & t )
	{
		return ( t > s );
	}


	/// @brief string >= Fstring
	friend
	inline
	bool
	operator >=( std::string const & s, Fstring const & t )
	{
		return ( t <= s );
	}


	/// @brief string > Fstring
	friend
	inline
	bool
	operator >( std::string const & s, Fstring const & t )
	{
		return ( t < s );
	}


	/// @brief Fstring <= cstring
	friend
	bool
	operator <=( Fstring const & s, c_cstring const t );


	/// @brief Fstring < cstring
	friend
	bool
	operator <( Fstring const & s, c_cstring const t );


	/// @brief Fstring >= cstring
	friend
	inline
	bool
	operator >=( Fstring const & s, c_cstring const t )
	{
		return !( s < t );
	}


	/// @brief Fstring > cstring
	friend
	inline
	bool
	operator >( Fstring const & s, c_cstring const t )
	{
		return !( s <= t );
	}


	/// @brief cstring <= Fstring
	friend
	inline
	bool
	operator <=( c_cstring const s, Fstring const & t )
	{
		return ( t >= s );
	}


	/// @brief cstring < Fstring
	friend
	inline
	bool
	operator <( c_cstring const s, Fstring const & t )
	{
		return ( t > s );
	}


	/// @brief cstring >= Fstring
	friend
	inline
	bool
	operator >=( c_cstring const s, Fstring const & t )
	{
		return ( t <= s );
	}


	/// @brief cstring > Fstring
	friend
	inline
	bool
	operator >( c_cstring const s, Fstring const & t )
	{
		return ( t < s );
	}


public: // Substring


	/// @brief Constant Substring: s( i, j )
	Fsubstring const
	operator ()( size_type const i, size_type const j ) const;


	/// @brief Substring: s( i, j )
	Fsubstring
	operator ()( size_type const i, size_type const j );


	/// @brief Constant Tail Substring: s( i )
	Fstring const
	operator ()( size_type const i ) const;


	/// @brief Tail Substring: s( i )
	Fsubstring
	operator ()( size_type const i );


	/// @brief Space-Free Head Constant Substring
	Fsubstring const
	head() const;


	/// @brief Space-Free Head Substring
	Fsubstring
	head();


	/// @brief Space Tail Substring
	Fsubstring
	tail();


	/// @brief Space Tail Constant Substring
	Fsubstring const
	tail() const;


public: // I/O


	/// @brief Stream Input
	friend
	std::istream &
	operator >>( std::istream & stream, Fstring & s );


	/// @brief Get from Stream
	friend
	std::istream &
	get( std::istream & stream, Fstring & s );


	/// @brief Get Line from Stream
	friend
	std::istream &
	getline( std::istream & stream, Fstring & s );


	/// @brief Read from Stream
	friend
	std::istream &
	read( std::istream & stream, Fstring & s );


	/// @brief Read Available Characters from Stream
	friend
	std::istream &
	readsome( std::istream & stream, Fstring & s );


	/// @brief Stream Output
	friend
	std::ostream &
	operator <<( std::ostream & stream, Fstring const & s );


private: // Data


	/// @brief Length
	size_type len_;

	/// @brief String
	char * str_;

	/// @brief cstring
	mutable char * c_str_;

	/// @brief Substring flag
	bool const sub_;


}; // Fstring


/// @brief Fstring + Fstring
Fstring
operator +( Fstring const & s, Fstring const & t );


/// @brief Fstring + string
std::string
operator +( Fstring const & s, std::string const & t );


/// @brief string + Fstring
#ifndef __clang__
std::string
operator +( std::string const & t, Fstring const & s );
#endif


/// @brief Fstring + cstring
Fstring
operator +( Fstring const & s, c_cstring const t );


/// @brief cstring + Fstring
Fstring
operator +( c_cstring const s, Fstring const & t );


/// @brief Fstring + char
Fstring
operator +( Fstring const & s, char const c );


/// @brief char + Fstring
Fstring
operator +( char const c, Fstring const & s );


/// @brief Fstring == Fstring
bool
operator ==( Fstring const & s, Fstring const & t );


/// @brief Fstring != Fstring
bool
operator !=( Fstring const & s, Fstring const & t );


/// @brief Fstring == string
bool
operator ==( Fstring const & s, std::string const & t );


/// @brief string == Fstring
bool
operator ==( std::string const & t, Fstring const & s );


/// @brief Fstring != string
bool
operator !=( Fstring const & s, std::string const & t );


/// @brief string != Fstring
bool
operator !=( std::string const & t, Fstring const & s );


/// @brief Fstring == cstring
bool
operator ==( Fstring const & s, c_cstring const t );


/// @brief cstring == Fstring
bool
operator ==( c_cstring const t, Fstring const & s );


/// @brief Fstring != cstring
bool
operator !=( Fstring const & s, c_cstring const t );


/// @brief cstring != Fstring
bool
operator !=( c_cstring const t, Fstring const & s );


/// @brief Fstring == char
bool
operator ==( Fstring const & s, char const c );


/// @brief char == Fstring
bool
operator ==( char const c, Fstring const & s );


/// @brief Fstring != char
bool
operator !=( Fstring const & s, char const c );


/// @brief char != Fstring
bool
operator !=( char const c, Fstring const & s );


/// @brief Fstring == Fstring Case-Insensitively?
bool
equali( Fstring const & s, Fstring const & t );


/// @brief Fstring == string Case-Insensitively?
bool
equali( Fstring const & s, std::string const & t );


/// @brief string == Fstring Case-Insensitively?
bool
equali( std::string const & s, Fstring const & t );


/// @brief Fstring == char Case-Insensitively?
bool
equali( Fstring const & s, char const c );


/// @brief char == Fstring Case-Insensitively?
bool
equali( char const c, Fstring const & s );


/// @brief Fstring == Fstring Case-Optionally?
bool
equal( Fstring const & s, Fstring const & t, bool const exact_case );


/// @brief Fstring == char Case-Optionally?
bool
equal( Fstring const & s, char const c, bool const exact_case );


/// @brief char == Fstring Case-Optionally?
bool
equal( char const c, Fstring const & s, bool const exact_case );


/// @brief Fstring <= Fstring
bool
operator <=( Fstring const & s, Fstring const & t );


/// @brief Fstring < Fstring
bool
operator <( Fstring const & s, Fstring const & t );


/// @brief Fstring >= Fstring
bool
operator >=( Fstring const & s, Fstring const & t );


/// @brief Fstring > Fstring
bool
operator >( Fstring const & s, Fstring const & t );


/// @brief Fstring <= string
bool
operator <=( Fstring const & s, std::string const & t );


/// @brief Fstring < string
bool
operator <( Fstring const & s, std::string const & t );


/// @brief Fstring >= string
bool
operator >=( Fstring const & s, std::string const & t );


/// @brief Fstring > string
bool
operator >( Fstring const & s, std::string const & t );


/// @brief string <= Fstring
bool
operator <=( std::string const & s, Fstring const & t );


/// @brief string < Fstring
bool
operator <( std::string const & s, Fstring const & t );


/// @brief string >= Fstring
bool
operator >=( std::string const & s, Fstring const & t );


/// @brief string > Fstring
bool
operator >( std::string const & s, Fstring const & t );


/// @brief Fstring <= cstring
bool
operator <=( Fstring const & s, c_cstring const t );


/// @brief Fstring < cstring
bool
operator <( Fstring const & s, c_cstring const t );


/// @brief Fstring >= cstring
bool
operator >=( Fstring const & s, c_cstring const t );


/// @brief Fstring > cstring
bool
operator >( Fstring const & s, c_cstring const t );


/// @brief cstring <= Fstring
bool
operator <=( c_cstring const s, Fstring const & t );


/// @brief cstring < Fstring
bool
operator <( c_cstring const s, Fstring const & t );


/// @brief cstring >= Fstring
bool
operator >=( c_cstring const s, Fstring const & t );


/// @brief cstring > Fstring
bool
operator >( c_cstring const s, Fstring const & t );


/// @brief Stream Input
std::istream &
operator >>( std::istream & stream, Fstring & s );


/// @brief Get from Stream
std::istream &
get( std::istream & stream, Fstring & s );


/// @brief Get Line from Stream
std::istream &
getline( std::istream & stream, Fstring & s );


/// @brief Read from Stream
std::istream &
read( std::istream & stream, Fstring & s );


/// @brief Read Available Characters from Stream
std::istream &
readsome( std::istream & stream, Fstring & s );


/// @brief Stream Output
std::ostream &
operator <<( std::ostream & stream, Fstring const & s );


// Fstring Member Function Explicit Specializations


	/// @brief Fstring is Readable as a char Supporting Stream Input?
	template<>
	inline
	bool
	Fstring::is_type< char >() const
	{
		return ( size() == 1 );
	}


	/// @brief Fstring is Readable as a string Supporting Stream Input?
	template<>
	inline
	bool
	Fstring::is_type< std::string >() const
	{
		return true;
	}


	/// @brief char of an Fstring
	template<>
	inline
	char
	Fstring::type_of< char >() const
	{
		return ( len_ == 1 ? str_[ 0 ] : char() ); // Check is_type first
	}


	/// @brief string of an Fstring
	template<>
	inline
	std::string
	Fstring::type_of< std::string >() const
	{
		return std::string( str_, len_ );
	}


/// @brief Fsubstring: Fixed-Length Fortran-Compatible Substring
///
/// @remarks
///  @li Subscripts run from 1 to the length
///  @li Space-padding is used in comparisons and assignments
///  @li Internal string rep is not null-terminated
///  @li Zero-length Fsubstrings are supported but cannot be indexed into (no valid indices)
///  @li Fsubstring not for explicit use in client code: Client code uses Fstring::operator () to get substrings
///  @li Pass s( i, j ).ref() to a non-const Fstring& argument
///  @li Don't return a substring of a local as an Fsubstring since its copy ctor uses ref semantics: Return as an Fstring to get a copy
class Fsubstring :
	public Fstring
{


private: // Types


	typedef  Fstring  Super;


	friend class Fstring;


public: // Creation


	/// @brief Copy Constructor
	inline
	Fsubstring( Fsubstring const & s ) :
		Fstring( s, 1, s.len_ )
	{}


	/// @brief Destructor
	inline
	virtual
	~Fsubstring()
	{}


private: // Creation


	/// @brief Fstring Range Constructor
	inline
	Fsubstring( Fstring const & s, size_type const i, size_type const j ) :
		Fstring( s, i, j )
	{}


	/// @brief Fstring Tail Constructor
	inline
	Fsubstring( Fstring const & s, size_type const i ) :
		Fstring( s, i )
	{}


public: // Assignment


	/// @brief Copy Assignment
	Fsubstring &
	operator =( Fsubstring const & s );


	/// @brief = Fstring
	Fsubstring &
	operator =( Fstring const & s );


	/// @brief = string
	Fsubstring &
	operator =( std::string const & s );


	/// @brief = cstring
	Fsubstring &
	operator =( c_cstring const s );


	/// @brief = char
	Fsubstring &
	operator =( char const c );


public: // Modifier


	/// @brief Reference to Fstring: Can Pass s( i, j ).ref() to an Fstring& Argument
	inline
	Fstring &
	ref()
	{
		return *this;
	}


}; // Fsubstring


// Fortran-Intrinsic-Compatible String Functions


namespace Fortran { // Control collision with Windows CHAR type and get Fstring, not char, output


/// @brief One-Character Fstring of a Given ASCII Integer Value
inline
Fstring
CHAR( int const i )
{
	return Fstring( static_cast< char >( i ) );
}


/// @brief One-Character Fstring of a Given ASCII Integer Value
inline
Fstring
ACHAR( int const i )
{
	return Fstring( static_cast< char >( i ) );
}


} // namespace Fortran


/// @brief Integer Value of a Given One-Character Fstring
inline
int
ICHAR( Fstring const & s )
{
	assert( s.length() == 1 );
	return static_cast< int >( s[ 1 ] );
}


/// @brief ASCII Integer Value of a Given One-Character Fstring
inline
int
IACHAR( Fstring const & s )
{
	assert( s.length() == 1 );
	return static_cast< int >( s[ 1 ] );
}


/// @brief First Index Position of a Substring in an Fstring
inline
Fstring::size_type
index( Fstring const & s, Fstring const & ss )
{
	return s.find( ss );
}


/// @brief First Index Position of a C-Substring in an Fstring
inline
Fstring::size_type
index( Fstring const & s, c_cstring const ss )
{
	return s.find( ss );
}


/// @brief First Index Position of a Character in an Fstring
inline
Fstring::size_type
index( Fstring const & s, char const c )
{
	return s.find( c );
}


/// @brief Length
inline
Fstring::size_type
len( Fstring const & s )
{
	return s.length();
}


/// @brief Length Space-Trimmed
inline
Fstring::size_type
len_trim( Fstring const & s )
{
	return s.len_trim();
}


/// @brief Space-Trimmed Copy
inline
Fstring
trimmed( Fstring const & s )
{
	return s.trimmed();
}


/// @brief ASCII Lexical >= Comparison
inline
bool
lge( Fstring const & s, Fstring const & t )
{
	return ( s >= t );
}


/// @brief ASCII Lexical < Comparison
inline
bool
lgt( Fstring const & s, Fstring const & t )
{
	return ( s > t );
}


/// @brief ASCII Lexical <= Comparison
inline
bool
lle( Fstring const & s, Fstring const & t )
{
	return ( s <= t );
}


/// @brief ASCII Lexical < Comparison
inline
bool
llt( Fstring const & s, Fstring const & t )
{
	return ( s < t );
}


// Fortran Migration Support String Functions


// Predicate


/// @brief Fstring is Blank?
inline
bool
is_blank( Fstring const & s )
{
	return s.is_blank();
}


/// @brief Fstring is Not Blank?
inline
bool
not_blank( Fstring const & s )
{
	return s.not_blank();
}


/// @brief Fstring Has Any Characters of a Set?
inline
bool
has_any_of( Fstring const & s, Fstring const & t )
{
	return s.has_any_of( t );
}


/// @brief Fstring Has Any Characters of a Set?
inline
bool
has_any_of( Fstring const & s, c_cstring const t )
{
	return s.has_any_of( t );
}


// Search


/// @brief Last Index Position of a Substring in an Fstring
inline
Fstring::size_type
last_index( Fstring const & s, Fstring const & ss )
{
	return s.find_last( ss );
}


/// @brief Last Index Position of a Substring in an Fstring
inline
Fstring::size_type
last_index( Fstring const & s, c_cstring const ss )
{
	return s.find_last( ss );
}


/// @brief Last Index Position of a Character in an Fstring
inline
Fstring::size_type
last_index( Fstring const & s, char const c )
{
	return s.find_last( c );
}


// Modifier


/// @brief Lowercase an Fstring
inline
Fstring &
lowercase( Fstring & s )
{
	return s.lowercase();
}


/// @brief Uppercase an Fstring
inline
Fstring &
uppercase( Fstring & s )
{
	return s.uppercase();
}


/// @brief Lowercase an Fstring
inline
void
str_dn( Fstring & s )
{
	s.lowercase();
}


/// @brief Uppercase an Fstring
inline
void
str_up( Fstring & s )
{
	s.uppercase();
}


/// @brief Lowercased Copy in an Output Fstring
inline
void
str_dncase( Fstring & s_out, Fstring const & s_in )
{
	s_out = s_in;
	s_out.lowercase();
}


/// @brief Uppercased Copy in an Output Fstring
inline
void
str_upcase( Fstring & s_out, Fstring const & s_in )
{
	s_out = s_in;
	s_out.uppercase();
}


// Generator


/// @brief Left-Justified Copy
inline
Fstring
ljust( Fstring const & s )
{
	return s.left_justified();
}


/// @brief Right-Justified Copy
inline
Fstring
rjust( Fstring const & s )
{
	return s.right_justified();
}


/// @brief Compressed Copy
inline
Fstring
compress( Fstring const & s )
{
	return s.compressed();
}


/// @brief Centered Copy
inline
Fstring
center( Fstring const & s )
{
	return s.centered();
}


/// @brief Lowercased Copy
inline
Fstring
lowercased( Fstring const & s )
{
	return s.lowercased();
}


/// @brief Uppercased Copy
inline
Fstring
uppercased( Fstring const & s )
{
	return s.uppercased();
}


/// @brief Lowercased Copy
inline
Fstring
dncase( Fstring const & s )
{
	return s.lowercased();
}


/// @brief Uppercased Copy
inline
Fstring
upcase( Fstring const & s )
{
	return s.uppercased();
}


// Conversion To Fstring


/// @brief Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
Fstring_of( T const & t )
{
	std::ostringstream t_stream;
	t_stream << std::uppercase << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief Fstring of a string Specialization
template<>
inline
Fstring
Fstring_of< std::string >( std::string const & t )
{
	return Fstring( t );
}


/// @brief Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
Fstring_of(
	T const & t,
	int const p // Precision
)
{
	std::ostringstream t_stream;
	t_stream << std::uppercase << std::setprecision( p ) << t;
	return t_stream.str();
}


/// @brief Left-Justified Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
left_Fstring_of(
	T const & t,
	int const w, // Minimum width
	char const f = ' ' // Fill character
)
{
	std::ostringstream t_stream;
	t_stream << std::left << std::uppercase
	 << std::setw( w ) << std::setfill( f ) << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief Right-Justified Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
right_Fstring_of(
	T const & t,
	int const w, // Minimum width
	char const f = ' ' // Fill character
)
{
	std::ostringstream t_stream;
	t_stream << std::right << std::uppercase
	 << std::setw( w ) << std::setfill( f ) << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief Leading-Zero Right-Justified Fstring of a Template Argument Type Supporting Stream Output
/// @note Negative numbers appear with the minus sign on the left of the filled zeros
template< typename T >
inline
Fstring
lead_zero_Fstring_of(
	T const & t,
	int const w // Minimum width
)
{
	std::ostringstream t_stream;
	t_stream << std::internal << std::uppercase
	 << std::setw( w ) << std::setfill( '0' ) << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief Right-Justified General Format Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
general_Fstring_of(
	T const & t,
	int const w = TypeTraits< T >::width(), // Minimum width
	int const p = TypeTraits< T >::precision() // Precision
)
{
	std::ostringstream t_stream;
	t_stream << std::right << std::uppercase << std::showpoint
	 << std::setw( w ) << std::setprecision( p ) << t;
	return t_stream.str();
}


/// @brief Right-Justified Fixed Format Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
fixed_Fstring_of(
	T const & t,
	int const w = TypeTraits< T >::width(), // Minimum width
	int const p = TypeTraits< T >::precision() // Precision
)
{
	std::ostringstream t_stream;
	t_stream << std::right << std::uppercase << std::fixed << std::showpoint
	 << std::setw( w ) << std::setprecision( p ) << t;
	return t_stream.str();
}


/// @brief Right-Justified Scientific Format Fstring of a Template Argument Type Supporting Stream Output
template< typename T >
inline
Fstring
scientific_Fstring_of(
	T const & t,
	int const w = TypeTraits< T >::width(), // Minimum width
	int const p = TypeTraits< T >::precision() // Precision
)
{
	std::ostringstream t_stream;
	t_stream << std::right << std::uppercase << std::scientific << std::showpoint
	 << std::setw( w ) << std::setprecision( p ) << t;
	return t_stream.str();
}


// Conversion From Fstring


/// @brief Fstring is Readable as a Type Supporting Stream Input?
template< typename T >
inline
bool
is_type( Fstring const & s )
{
	if ( s.is_blank() ) { // Don't accept blank Fstring
		return false;
	} else { // Try to read the Fstring as a T
		std::istringstream t_stream( s.trimmed_whitespace() );
		T t;
		t_stream >> t;
		return ( ( t_stream ) && ( t_stream.eof() ) );
	}
}


/// @brief Fstring is Readable as a string Supporting Stream Input?
template<>
inline
bool
is_type< std::string >( Fstring const & )
{
	return true;
}


/// @brief Fstring is Readable as a char Supporting Stream Input?
template<>
inline
bool
is_type< char >( Fstring const & s )
{
	return ( s.size() == 1 );
}


/// @brief Fstring is Readable as a bool?
inline
bool
is_bool( Fstring const & s )
{
	return is_type< bool >( s );
}


/// @brief Fstring is Readable as a short int?
inline
bool
is_short( Fstring const & s )
{
	return is_type< short int >( s );
}


/// @brief Fstring is Readable as an int?
inline
bool
is_int( Fstring const & s )
{
	return is_type< int >( s );
}


/// @brief Fstring is Readable as a long int?
inline
bool
is_long( Fstring const & s )
{
	return is_type< long int >( s );
}


/// @brief Fstring is Readable as a unsigned short int?
inline
bool
is_ushort( Fstring const & s )
{
	return is_type< unsigned short int >( s );
}


/// @brief Fstring is Readable as an unsigned int?
inline
bool
is_uint( Fstring const & s )
{
	return is_type< unsigned int >( s );
}


/// @brief Fstring is Readable as a unsigned long int?
inline
bool
is_ulong( Fstring const & s )
{
	return is_type< unsigned long int >( s );
}


/// @brief Fstring is Readable as a float?
inline
bool
is_float( Fstring const & s )
{
	return is_type< float >( s );
}


/// @brief Fstring is Readable as a double?
inline
bool
is_double( Fstring const & s )
{
	return is_type< double >( s );
}


/// @brief Fstring is Readable as a long double?
inline
bool
is_longdouble( Fstring const & s )
{
	return is_type< long double >( s );
}


/// @brief Fstring is Readable as a char?
inline
bool
is_char( Fstring const & s )
{
	return is_type< char >( s );
}


/// @brief Fstring is Readable as a string?
inline
bool
is_string( Fstring const & )
{
	return true;
}


/// @brief Type of an Fstring for Type Supporting Stream Input
template< typename T >
inline
T
type_of( Fstring const & s )
{
	std::istringstream t_stream( s.trimmed_whitespace() );
	T t;
	t_stream >> t;
	return ( ( t_stream ) && ( t_stream.eof() ) ? t : T() ); // Check is_type first
}


/// @brief string of an Fstring
template<>
inline
std::string
type_of< std::string >( Fstring const & s )
{
	return std::string( s );
}


/// @brief char of an Fstring
template<>
inline
char
type_of< char >( Fstring const & s )
{
	return ( s.size() == 1 ? s[ 0 ] : char() ); // Check is_type first
}


/// @brief short int of an Fstring
inline
short int
short_of( Fstring const & s )
{
	return type_of< short int >( s );
}


/// @brief int of an Fstring
inline
int
int_of( Fstring const & s )
{
	return type_of< int >( s );
}


/// @brief long int of an Fstring
inline
long int
long_of( Fstring const & s )
{
	return type_of< long int >( s );
}


/// @brief unsigned short int of an Fstring
inline
unsigned short int
ushort_of( Fstring const & s )
{
	return type_of< unsigned short int >( s );
}


/// @brief unsigned int of an Fstring
inline
unsigned int
uint_of( Fstring const & s )
{
	return type_of< unsigned int >( s );
}


/// @brief unsigned long int of an Fstring
inline
unsigned long int
ulong_of( Fstring const & s )
{
	return type_of< unsigned long int >( s );
}


/// @brief float of an Fstring
inline
float
float_of( Fstring const & s )
{
	return type_of< float >( s );
}


/// @brief double of an Fstring
inline
double
double_of( Fstring const & s )
{
	return type_of< double >( s );
}


/// @brief long double of an Fstring
inline
long double
longdouble_of( Fstring const & s )
{
	return type_of< long double >( s );
}


/// @brief char of an Fstring
inline
char
char_of( Fstring const & s )
{
	return type_of< char >( s );
}


/// @brief string of an Fstring
inline
std::string
string_of( Fstring const & s )
{
	return std::string( s );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_Fstring_HH
