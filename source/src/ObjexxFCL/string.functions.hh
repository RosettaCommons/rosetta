#ifndef INCLUDED_ObjexxFCL_string_functions_HH
#define INCLUDED_ObjexxFCL_string_functions_HH


// String Functions
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
#include <ObjexxFCL/TypeTraits.hh>

// C++ Headers
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>


namespace ObjexxFCL {


// Types
typedef  char       *    cstring;
typedef  char const *  c_cstring;


// Predicate


/// @brief string == string Case-Insensitively?
bool
equali( std::string const & s, std::string const & t );


/// @brief string == cstring Case-Insensitively?
bool
equali( std::string const & s, c_cstring const t );


/// @brief cstring == string Case-Insensitively?
bool
equali( c_cstring const s, std::string const & t );


/// @brief string == string Case-Optionally?
inline
bool
equal( std::string const & s, std::string const & t, bool const exact_case = true )
{
	if ( exact_case ) {
		return ( s == t );
	} else {
		return equali( s, t );
	}
}


/// @brief string is Blank?
inline
bool
is_blank( std::string const & s )
{
	if ( s.empty() ) {
		return true;
	} else {
		return ( s.find_first_not_of( ' ' ) == std::string::npos );
	}
}


/// @brief string is Not Blank?
inline
bool
not_blank( std::string const & s )
{
	return ( ! is_blank( s ) );
}


/// @brief string is Whitespace?
inline
bool
is_whitespace( std::string const & s )
{
	if ( s.empty() ) {
		return true;
	} else {
		return ( s.find_last_not_of( " \t\000" ) == std::string::npos );
	}
}


/// @brief string is Not Whitespace?
inline
bool
not_whitespace( std::string const & s )
{
	return ( ! is_whitespace( s ) );
}


/// @brief string has a string?
inline
bool
has( std::string const & s, std::string const & t )
{
	return ( s.find( t ) != std::string::npos );
}


/// @brief string has a cstring?
inline
bool
has( std::string const & s, c_cstring const t )
{
	return ( s.find( t ) != std::string::npos );
}


/// @brief string has a Character?
inline
bool
has( std::string const & s, char const c )
{
	return ( s.find( c ) != std::string::npos );
}


/// @brief string has Any Character of a string?
inline
bool
has_any_of( std::string const & s, std::string const & t )
{
	return ( s.find_first_of( t ) != std::string::npos );
}


/// @brief string has a Character?
inline
bool
has_any_of( std::string const & s, char const c )
{
	return ( s.find( c ) != std::string::npos );
}


/// @brief Has a Prefix Case-Optionally?
bool
has_prefix( std::string const & s, std::string const & pre, bool const exact_case = true );


/// @brief Has a Suffix Case-Optionally?
bool
has_suffix( std::string const & s, std::string const & suf, bool const exact_case = true );


// Inspector


/// @brief Length Space-Trimmed
inline
std::string::size_type
len_trim( std::string const & s )
{
	return s.find_last_not_of( ' ' ) + 1; // Works if npos returned: npos + 1 == 0
}


/// @brief Length Whitespace-Trimmed
inline
std::string::size_type
len_trim_whitespace( std::string const & s )
{
	return s.find_last_not_of( " \t\000" ) + 1; // Works if npos returned: npos + 1 == 0
}


// Modifier


/// @brief Lowercase a string
std::string &
lowercase( std::string & s );


/// @brief Uppercase a string
std::string &
uppercase( std::string & s );


/// @brief Left Justify a string
std::string &
left_justify( std::string & s );


/// @brief Right Justify a string
std::string &
right_justify( std::string & s );


/// @brief Trim Trailing Space from a string
std::string &
trim( std::string & s );


/// @brief Trim Trailing Whitespace from a string
std::string &
trim_whitespace( std::string & s );


/// @brief Strip Specified Characters from a string's Tails
std::string &
strip( std::string & s, std::string const & chars );


/// @brief Strip Specified Characters from a string's Left Tail
std::string &
lstrip( std::string & s, std::string const & chars );


/// @brief Strip Specified Characters from a string's Right Tail
std::string &
rstrip( std::string & s, std::string const & chars );


/// @brief Strip Space from a string's Tails
std::string &
strip( std::string & s );


/// @brief Strip Space from a string's Left Tail
std::string &
lstrip( std::string & s );


/// @brief Strip Space from a string's Right Tail
std::string &
rstrip( std::string & s );


/// @brief Strip Whitespace from a string's Tails
std::string &
strip_whitespace( std::string & s );


/// @brief Strip Whitespace from a string's Left Tail
std::string &
lstrip_whitespace( std::string & s );


/// @brief Strip Whitespace from a string's Right Tail
std::string &
rstrip_whitespace( std::string & s );


/// @brief Pad a string to a Specified Length
std::string &
pad( std::string & s, std::string::size_type len );


/// @brief Left-Pad a string to a Specified Length
std::string &
lpad( std::string & s, std::string::size_type len );


/// @brief Right-Pad a string to a Specified Length
std::string &
rpad( std::string & s, std::string::size_type len );


/// @brief Size a string to a Specified Length
std::string &
size( std::string & s, std::string::size_type len );


/// @brief Center a string wrt its Whitespace
std::string &
center( std::string & s );


/// @brief Center a string with a Specified Length
std::string &
center( std::string & s, std::string::size_type len );


/// @brief Remove Repeat Characters from a Possibly Unsorted string Preserving Order
std::string &
unique( std::string & s );


/// @brief Overlay a string With Another string, Expanding Size as Needed
std::string &
overlay( std::string & s, std::string const & t, std::string::size_type pos = 0 );


// Generator


/// @brief Blank string of Specified Length
inline
std::string
blank( std::string::size_type len )
{
	return std::string( len, ' ' );
}


/// @brief Lowercased Copy of a string
std::string
lowercased( std::string const & s );


/// @brief Uppercased Copy of a string
std::string
uppercased( std::string const & s );


/// @brief Left-Justified Copy of a string
std::string
left_justified( std::string const & s );


/// @brief Right-Justified Copy of a string
std::string
right_justified( std::string const & s );


/// @brief Trailing Space Trimmed Copy of a string
std::string
trimmed( std::string const & s );


/// @brief Trailing Whitespace Trimmed Copy of a string
std::string
trimmed_whitespace( std::string const & s );


/// @brief Specified Characters Stripped from a string's Tails Copy of a string
std::string
stripped( std::string const & s, std::string const & chars );


/// @brief Specified Characters Stripped from a string's Left Tail Copy of a string
std::string
lstripped( std::string const & s, std::string const & chars );


/// @brief Specified Characters Stripped from a string's Right Tail Copy of a string
std::string
rstripped( std::string const & s, std::string const & chars );


/// @brief Space Stripped from a string's Tails Copy of a string
std::string
stripped( std::string const & s );


/// @brief Space Stripped from a string's Left Tail Copy of a string
std::string
lstripped( std::string const & s );


/// @brief Space Stripped from a string's Right Tail Copy of a string
std::string
rstripped( std::string const & s );


/// @brief Whitespace Stripped from a string's Tails Copy of a string
std::string
stripped_whitespace( std::string const & s );


/// @brief Whitespace Stripped from a string's Left Tail Copy of a string
std::string
lstripped_whitespace( std::string const & s );


/// @brief Whitespace Stripped from a string's Right Tail Copy of a string
std::string
rstripped_whitespace( std::string const & s );


/// @brief Padded to a Specified Length Copy of a string
std::string
padded( std::string const & s, std::string::size_type len );


/// @brief Left-Padded to a Specified Length Copy of a string
std::string
lpadded( std::string const & s, std::string::size_type len );


/// @brief Right-Padded to a Specified Length Copy of a string
std::string
rpadded( std::string const & s, std::string::size_type len );


/// @brief Sized to a Specified Length Copy of a string
std::string
sized( std::string const & s, std::string::size_type len );


/// @brief Centered wrt Whitespace Copy of a string
std::string
centered( std::string const & s );


/// @brief Centered in a string of Specified Length Copy of a string
std::string
centered( std::string const & s, std::string::size_type len );


/// @brief Removed Repeat Characters from a Possibly Unsorted string Preserving Order Copy of a string
std::string
uniqued( std::string const & s );


/// @brief Space-Free Head Copy of a string
std::string
head( std::string const & s );


// Conversion To std::string


/// @brief string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
string_of( T const & t )
{
	std::ostringstream t_stream;
	t_stream << std::uppercase << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
string_of(
	T const & t,
	int const p // Precision
)
{
	std::ostringstream t_stream;
	t_stream << std::uppercase << std::setprecision( p ) << t;
	return t_stream.str();
}


/// @brief Left-Justified string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
left_string_of(
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


/// @brief Right-Justified string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
right_string_of(
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


/// @brief Leading-Zero Right-Justified string of a Template Argument Type Supporting Stream Output
/// @note Negative numbers appear with the minus sign on the left of the filled zeros
template< typename T >
inline
std::string
lead_zero_string_of(
	T const & t,
	int const w // Minimum width
)
{
	std::ostringstream t_stream;
	t_stream << std::internal << std::uppercase
	 << std::setw( w ) << std::setfill( '0' ) << std::setprecision( TypeTraits< T >::precision() ) << t;
	return t_stream.str();
}


/// @brief Right-Justified General Format string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
general_string_of(
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


/// @brief Right-Justified Fixed Format string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
fixed_string_of(
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


/// @brief Right-Justified Scientific Format string of a Template Argument Type Supporting Stream Output
template< typename T >
inline
std::string
scientific_string_of(
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


// Conversion From std::string


/// @brief string is Readable as a Type Supporting Stream Input?
template< typename T >
inline
bool
is_type( std::string const & s )
{
	if ( is_whitespace( s ) ) { // Don't accept empty or whitespace string
		return false;
	} else { // Try to read the string as a T
		std::istringstream t_stream( trimmed_whitespace( s ) );
		T t;
		t_stream >> t;
		return ( ( t_stream ) && ( t_stream.eof() ) );
	}
}


/// @brief string is Readable as a char Supporting Stream Input?
template<>
inline
bool
is_type< char >( std::string const & s )
{
	return ( s.length() == 1 );
}


/// @brief string is Readable as a bool?
inline
bool
is_bool( std::string const & s )
{
	return is_type< bool >( s );
}


/// @brief string is Readable as a short int?
inline
bool
is_short( std::string const & s )
{
	return is_type< short int >( s );
}


/// @brief string is Readable as an int?
inline
bool
is_int( std::string const & s )
{
	return is_type< int >( s );
}

/// @brief string is Readable as ints? [e.g., "5" or "5-8"]
bool
is_ints( std::string const & s );

/// @brief string is Readable as a long int?
inline
bool
is_long( std::string const & s )
{
	return is_type< long int >( s );
}


/// @brief string is Readable as a unsigned short int?
inline
bool
is_ushort( std::string const & s )
{
	return is_type< unsigned short int >( s );
}


/// @brief string is Readable as an unsigned int?
inline
bool
is_uint( std::string const & s )
{
	return is_type< unsigned int >( s );
}


/// @brief string is Readable as a unsigned long int?
inline
bool
is_ulong( std::string const & s )
{
	return is_type< unsigned long int >( s );
}


/// @brief string is Readable as a float?
inline
bool
is_float( std::string const & s )
{
	return is_type< float >( s );
}


/// @brief string is Readable as a double?
inline
bool
is_double( std::string const & s )
{
	return is_type< double >( s );
}


/// @brief string is Readable as a long double?
inline
bool
is_longdouble( std::string const & s )
{
	return is_type< long double >( s );
}


/// @brief string is Readable as a char?
inline
bool
is_char( std::string const & s )
{
	return is_type< char >( s );
}


/// @brief Type of a string for Type Supporting Stream Input
template< typename T >
inline
T
type_of( std::string const & s )
{
	std::istringstream t_stream( trimmed_whitespace( s ) );
	T t;
	t_stream >> t;
	return ( ( t_stream ) && ( t_stream.eof() ) ? t : T() ); // Check is_type first
}


/// @brief char of a string
template<>
inline
char
type_of< char >( std::string const & s )
{
	return ( s.length() == 1 ? s[ 0 ] : char() ); // Check is_type first
}


/// @brief short int of a string
inline
short int
short_of( std::string const & s )
{
	return type_of< short int >( s );
}


/// @brief int of a string
inline
int
int_of( std::string const & s )
{
	return type_of< int >( s );
}

/// @brief ints of a string (e.g., allowing "5-8" to represent "5 6 7 8")
std::vector< int >
ints_of( std::string const & s );


/// @brief ints of a string (e.g., allowing "5-8" to represent "5 6 7 8")
std::vector< int >
ints_of( std::string const & s, bool & string_is_ok );


/// @brief long int of a string
inline
long int
long_of( std::string const & s )
{
	return type_of< long int >( s );
}


/// @brief unsigned short int of a string
inline
unsigned short int
ushort_of( std::string const & s )
{
	return type_of< unsigned short int >( s );
}


/// @brief unsigned int of a string
inline
unsigned int
uint_of( std::string const & s )
{
	return type_of< unsigned int >( s );
}


/// @brief unsigned long int of a string
inline
unsigned long int
ulong_of( std::string const & s )
{
	return type_of< unsigned long int >( s );
}


/// @brief float of a string
inline
float
float_of( std::string const & s )
{
	return type_of< float >( s );
}


/// @brief double of a string
inline
double
double_of( std::string const & s )
{
	return type_of< double >( s );
}


/// @brief long double of a string
inline
long double
longdouble_of( std::string const & s )
{
	return type_of< long double >( s );
}


/// @brief char of a string
inline
char
char_of( std::string const & s )
{
	return type_of< char >( s );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_string_functions_HH
