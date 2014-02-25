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
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <algorithm>
#include <cctype>
#include <cstring>


namespace ObjexxFCL {


// Using
using std::string;


// Constants
char const SPACE( ' ' );
std::string const WHITESPACE( " \t\0", 3 );


// Predicate


/// @brief char == char Case-Insensitively (non-inline for use by equali below)?
bool
char_equali( char const c, char const d )
{
	return ( std::tolower( c ) == std::tolower( d ) );
}


/// @brief string == string Case-Insensitively?
bool
equali( std::string const & s, std::string const & t )
{
	if ( s.length() != t.length() ) {
		return false;
	} else {
		return std::equal( s.begin(), s.end(), t.begin(), char_equali );
	}
}


/// @brief string == cstring Case-Insensitively?
bool
equali( std::string const & s, c_cstring const t )
{
	if ( s.length() != std::strlen( t ) ) {
		return false;
	} else {
		return std::equal( s.begin(), s.end(), t, char_equali );
	}
}


/// @brief cstring == string Case-Insensitively?
bool
equali( c_cstring const s, std::string const & t )
{
	if ( std::strlen( s ) != t.length() ) {
		return false;
	} else {
		return std::equal( t.begin(), t.end(), s, char_equali );
	}
}


/// @brief Has a Prefix Case-Optionally?
bool
has_prefix( std::string const & s, std::string const & pre, bool const exact_case )
{
	string::size_type const pre_len( pre.length() );
	if ( pre_len == 0 ) {
		return false;
	} else if ( s.length() < pre_len ) {
		return false;
	} else if ( exact_case ) {
		return ( s.find( pre ) == 0 );
	} else {
		return ( lowercased( s ).find( lowercased( pre ) ) == 0 );
	}
}


/// @brief Has a Suffix Case-Optionally?
bool
has_suffix( std::string const & s, std::string const & suf, bool const exact_case )
{
	string::size_type const suf_len( suf.length() );
	if ( suf_len == 0 ) {
		return false;
	} else {
		string::size_type const s_len( s.length() );
		if ( s_len < suf_len ) {
			return false;
		} else if ( exact_case ) {
			return ( s.rfind( suf ) == s_len - suf_len );
		} else {
			return ( lowercased( s ).rfind( lowercased( suf ) ) == s_len - suf_len );
		}
	}
}


// Modifier


/// @brief Lowercase a string
std::string &
lowercase( std::string & s )
{
	string::size_type const s_len( s.length() );
	for ( string::size_type i = 0; i < s_len; ++i ) {
		s[ i ] = std::tolower( s[ i ] );
	}
	return s;
}


/// @brief Uppercase a string
std::string &
uppercase( std::string & s )
{
	string::size_type const s_len( s.length() );
	for ( string::size_type i = 0; i < s_len; ++i ) {
		s[ i ] = std::toupper( s[ i ] );
	}
	return s;
}


/// @brief Left Justify a string
std::string &
left_justify( std::string & s )
{
	string::size_type const off( s.find_first_not_of( SPACE ) );
	if ( ( off > 0 ) && ( off != string::npos ) ) {
		s.erase( 0, off ).append( off, SPACE );
	}
	return s;
}


/// @brief Right Justify a string
std::string &
right_justify( std::string & s )
{
	string::size_type const s_len_trim( len_trim( s ) );
	string::size_type const off( s.length() - s_len_trim );
	if ( off > 0 ) {
		s.erase( s_len_trim ).insert( 0, off, SPACE );
	}
	return s;
}


/// @brief Trim Trailing Space from a string
std::string &
trim( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ie == string::npos ) { // Blank string: return empty string
			s.clear();
		} else if ( ie + 1 < s.length() ) { // Trim tail
			s.erase( ie + 1 );
		}
	}
	return s;
}


/// @brief Trim Trailing Whitespace from a string
std::string &
trim_whitespace( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ie == string::npos ) { // Blank string: return empty string
			s.clear();
		} else if ( ie + 1 < s.length() ) { // Trim tail
			s.erase( ie + 1 );
		}
	}
	return s;
}


/// @brief Strip Specified Characters from a string's Tails
std::string &
strip( std::string & s, std::string const & chars )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( chars ) );
		string::size_type const ie( s.find_last_not_of( chars ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is from chars
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
			if ( ib > 0 ) s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Specified Characters from a string's Left Tail
std::string &
lstrip( std::string & s, std::string const & chars )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( chars ) );
		if ( ib == string::npos ) { // All of string is from chars
			s.clear();
		} else if ( ib > 0 ) {
			s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Specified Characters from a string's Right Tail
std::string &
rstrip( std::string & s, std::string const & chars )
{
	if ( ! s.empty() ) {
		string::size_type const ie( s.find_last_not_of( chars ) );
		if ( ie == string::npos ) { // All of string is from chars
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
		}
	}
	return s;
}


/// @brief Strip Space from a string's Tails
std::string &
strip( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( SPACE ) );
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is SPACE
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
			if ( ib > 0 ) s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Space from a string's Left Tail
std::string &
lstrip( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( SPACE ) );
		if ( ib == string::npos ) { // All of string is SPACE
			s.clear();
		} else if ( ib > 0 ) {
			s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Space from a string's Right Tail
std::string &
rstrip( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ie == string::npos ) { // All of string is SPACE
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
		}
	}
	return s;
}


/// @brief Strip Whitespace from a string's Tails
std::string &
strip_whitespace( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( WHITESPACE ) );
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is from WHITESPACE
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
			if ( ib > 0 ) s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Whitespace from a string's Left Tail
std::string &
lstrip_whitespace( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ib( s.find_first_not_of( WHITESPACE ) );
		if ( ib == string::npos ) { // All of string is from WHITESPACE
			s.clear();
		} else if ( ib > 0 ) {
			s.erase( 0, ib );
		}
	}
	return s;
}


/// @brief Strip Whitespace from a string's Right Tail
std::string &
rstrip_whitespace( std::string & s )
{
	if ( ! s.empty() ) {
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ie == string::npos ) { // All of string is from WHITESPACE
			s.clear();
		} else {
			if ( ie < s.length() - 1 ) s.erase( ie + 1 );
		}
	}
	return s;
}


/// @brief Pad a string to a Specified Length
std::string &
pad( std::string & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Pad
		s.append( len - s_len, SPACE );
	}
	return s;
}


/// @brief Left-Pad a string to a Specified Length
std::string &
lpad( std::string & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Left-pad
		s.insert( static_cast< string::size_type >( 0 ), len - s_len, SPACE );
	}
	return s;
}


/// @brief Right-Pad a string to a Specified Length
std::string &
rpad( std::string & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Pad
		s.append( len - s_len, SPACE );
	}
	return s;
}


/// @brief Size a string to a Specified Length
std::string &
size( std::string & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Pad
		s.append( len - s_len, SPACE );
	} else if ( s_len > len ) { // Truncate
		s.erase( len );
	}
	return s;
}


/// @brief Center a string wrt its Whitespace
std::string &
center( std::string & s )
{
	string::size_type const s_len( s.length() );
	s = centered( strip_whitespace( s ), s_len );
	return s;
}


/// @brief Center a string with a Specified Length
std::string &
center( std::string & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Pad
		string::size_type const off( ( len - s_len ) / 2 );
		s = string( off, SPACE ).append( s ).append( string( len - s_len - off, SPACE ) );
	} else if ( s_len > len  ) { // Truncate
		s.erase( len );
	}
	return s;
}


/// @brief Remove Repeat Characters from a Possibly Unsorted string Preserving Order
std::string &
unique( std::string & s )
{
	string u;
	string::size_type const s_len( s.length() );
	for ( string::size_type i = 0; i < s_len; ++i ) {
		if ( u.find( s[ i ] ) == string::npos ) {
			u.push_back( s[ i ] );
		}
	}
	s.swap( u );
	return s;
}


/// @brief Overlay a string With Another string, Expanding Size as Needed
std::string &
overlay( std::string & s, std::string const & t, std::string::size_type const pos )
{
	std::string::size_type const t_len( t.length() );
	std::string::size_type const l_len( pos + t_len ); // Lower bound on new string length
	if ( l_len > s.length() ) s.resize( l_len, ' ' ); // Expand
	s.replace( pos, t_len, t ); // Overlay the string
	return s;
}


// Generator


/// @brief Lowercased Copy of a string
std::string
lowercased( std::string const & s )
{
	string t( s );
	string::size_type const t_len( t.length() );
	for ( string::size_type i = 0; i < t_len; ++i ) {
		t[ i ] = std::tolower( t[ i ] );
	}
	return t;
}


/// @brief Uppercased Copy of a string
std::string
uppercased( std::string const & s )
{
	string t( s );
	string::size_type const t_len( t.length() );
	for ( string::size_type i = 0; i < t_len; ++i ) {
		t[ i ] = std::toupper( t[ i ] );
	}
	return t;
}


/// @brief Left-Justified Copy of a string
std::string
left_justified( std::string const & s )
{
	string::size_type const off( s.find_first_not_of( SPACE ) );
	if ( ( off > 0 ) && ( off != string::npos ) ) {
		return s.substr( off ).append( off, SPACE );
	} else {
		return s;
	}
}


/// @brief Right-Justified Copy of a string
std::string
right_justified( std::string const & s )
{
	string::size_type const s_len_trim( len_trim( s ) );
	string::size_type const off( s.length() - s_len_trim );
	if ( off > 0 ) {
		return string( off, SPACE ).append( s.substr( 0, s_len_trim ) );
	} else {
		return s;
	}
}


/// @brief Trailing Space Trimmed Copy of a string
std::string
trimmed( std::string const & s )
{
	if ( s.empty() ) { // Empty string
		return s;
	} else {
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ie == string::npos ) { // Blank string: return empty string
			return string();
		} else if ( ie < s.length() - 1 ) { // Trimmed
			return s.substr( 0, ie + 1 );
		} else { // Unchanged
			return s;
		}
	}
}


/// @brief Trailing Whitespace Trimmed Copy of a string
std::string
trimmed_whitespace( std::string const & s )
{
	if ( s.empty() ) { // Empty string
		return s;
	} else {
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ie == string::npos ) { // Blank string: return empty string
			return string();
		} else if ( ie < s.length() - 1 ) { // Trimmed
			return s.substr( 0, ie + 1 );
		} else { // Unchanged
			return s;
		}
	}
}


/// @brief Specified Characters Stripped from a string's Tails Copy of a string
std::string
stripped( std::string const & s, std::string const & chars )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( chars ) );
		string::size_type const ie( s.find_last_not_of( chars ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is from chars
			return string(); // Return empty string
		} else {
			return s.substr( ib, ie - ib + 1 );
		}
	}
}


/// @brief Specified Characters Stripped from a string's Left Tail Copy of a string
std::string
lstripped( std::string const & s, std::string const & chars )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( chars ) );
		if ( ib == string::npos ) { // All of string is from chars
			return string(); // Return empty string
		} else if ( ib > 0 ) {
			return s.substr( ib );
		} else {
			return s;
		}
	}
}


/// @brief Specified Characters Stripped from a string's Right Tail Copy of a string
std::string
rstripped( std::string const & s, std::string const & chars )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ie( s.find_last_not_of( chars ) );
		if ( ie == string::npos ) { // All of string is from chars
			return string(); // Return empty string
		} else {
			return s.substr( 0, ie + 1 );
		}
	}
}


/// @brief Space Stripped from a string's Tails Copy of a string
std::string
stripped( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( SPACE ) );
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is SPACE
			return string(); // Return empty string
		} else {
			return s.substr( ib, ie - ib + 1 );
		}
	}
}


/// @brief Space Stripped from a string's Left Tail Copy of a string
std::string
lstripped( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( SPACE ) );
		if ( ib == string::npos ) { // All of string is SPACE
			return string(); // Return empty string
		} else if ( ib > 0 ) {
			return s.substr( ib );
		} else {
			return s;
		}
	}
}


/// @brief Space Stripped from a string's Right Tail Copy of a string
std::string
rstripped( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ie( s.find_last_not_of( SPACE ) );
		if ( ie == string::npos ) { // All of string is SPACE
			return string(); // Return empty string
		} else {
			return s.substr( 0, ie + 1 );
		}
	}
}


/// @brief Whitespace Stripped from a string's Tails Copy of a string
std::string
stripped_whitespace( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( WHITESPACE ) );
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ( ib == string::npos ) || ( ie == string::npos ) ) { // All of string is from WHITESPACE
			return string(); // Return empty string
		} else {
			return s.substr( ib, ie - ib + 1 );
		}
	}
}


/// @brief Whitespace Stripped from a string's Left Tail Copy of a string
std::string
lstripped_whitespace( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ib( s.find_first_not_of( WHITESPACE ) );
		if ( ib == string::npos ) { // All of string is from WHITESPACE
			return string(); // Return empty string
		} else if ( ib > 0 ) {
			return s.substr( ib );
		} else {
			return s;
		}
	}
}


/// @brief Whitespace Stripped from a string's Right Tail Copy of a string
std::string
rstripped_whitespace( std::string const & s )
{
	if ( s.empty() ) {
		return s;
	} else {
		string::size_type const ie( s.find_last_not_of( WHITESPACE ) );
		if ( ie == string::npos ) { // All of string is from WHITESPACE
			return string(); // Return empty string
		} else {
			return s.substr( 0, ie + 1 );
		}
	}
}


/// @brief Padded to a Specified Length Copy of a string
std::string
padded( std::string const & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Padded
		return s + string( len - s_len, SPACE );
	} else { // Unchanged
		return s;
	}
}


/// @brief Left-Padded to a Specified Length Copy of a string
std::string
lpadded( std::string const & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Left-padded
		return string( len - s_len, SPACE ).append( s );
	} else { // Unchanged
		return s;
	}
}


/// @brief Right-Padded to a Specified Length Copy of a string
std::string
rpadded( std::string const & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Padded
		return s + string( len - s_len, SPACE );
	} else { // Unchanged
		return s;
	}
}


/// @brief Sized to a Specified Length Copy of a string
std::string
sized( std::string const & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Padded
		return s + string( len - s_len, SPACE );
	} else if ( s_len == len  ) { // Unchanged
		return s;
	} else { // Truncated
		return s.substr( 0, len );
	}
}


/// @brief Centered wrt Whitespace Copy of a string
std::string
centered( std::string const & s )
{
	return centered( stripped_whitespace( s ), s.length() );
}


/// @brief Centered in a string of Specified Length Copy of a string
std::string
centered( std::string const & s, std::string::size_type const len )
{
	string::size_type const s_len( s.length() );
	if ( s_len < len ) { // Padded
		string::size_type const off( ( len - s_len ) / 2 );
		return string( off, SPACE ).append( s ).append( string( len - s_len - off, SPACE ) );
	} else if ( s_len == len  ) { // Unchanged
		return s;
	} else { // Truncated
		return s.substr( 0, len );
	}
}


/// @brief Removed Repeat Characters from a Possibly Unsorted string Preserving Order Copy of a string
std::string
uniqued( std::string const & s )
{
	string u;
	string::size_type const s_len( s.length() );
	for ( string::size_type i = 0; i < s_len; ++i ) {
		if ( u.find( s[ i ] ) == string::npos ) {
			u.push_back( s[ i ] );
		}
	}
	return u;
}


/// @brief Space-Free Head Copy of a string
std::string
head( std::string const & s )
{
	if ( s.empty() ) { // Empty string
		return s;
	} else {
		string::size_type const ie( s.find( SPACE ) );
		if ( ie == string::npos ) { // Space-free string
			return s;
		} else {
			return s.substr( 0, ie );
		}
	}
}


/// @brief  ints of a string (e.g., allowing "5-8" to represent "5 6 7 8").
bool
is_ints( std::string const & s ){
  bool string_is_ok( false );
  ints_of( s, string_is_ok );
  return string_is_ok;
}

/// @brief  ints of a string (e.g., allowing "5-8" to represent "5 6 7 8").
std::vector< int >
ints_of( std::string const & s ){
  bool string_is_ok( false );
  return ints_of( s, string_is_ok );
}

/// @brief  ints of a string (e.g., allowing "5-8" to represent "5 6 7 8").
std::vector< int >
ints_of( std::string const & s, bool & string_is_ok ){

  std::vector< int > vals;
  string_is_ok = false;

  size_t found_dash = s.find( "-" );
  if ( found_dash == 0 ) found_dash = s.substr(1).find("-") + 1;
  if ( found_dash == std::string::npos || found_dash == 0 ){

    string_is_ok = is_int( s );
    if ( string_is_ok ) vals.push_back( int_of( s ) );

  } else {
    std::string const start_val_string = s.substr(0,found_dash  );
    std::string const  end_val_string  = s.substr(found_dash+1, s.size() );
    string_is_ok = is_int( start_val_string ) && is_int( end_val_string); // currently cannot process
    if ( string_is_ok ){
      int start_val = int_of( start_val_string );
      int end_val = int_of( end_val_string );
      for ( int n = start_val; n <= end_val; n++ ) vals.push_back( n );
    }
  }

  return vals;
}


} // namespace ObjexxFCL
