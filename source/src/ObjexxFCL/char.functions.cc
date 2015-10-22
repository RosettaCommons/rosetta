// Character Functions
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
#include <ObjexxFCL/char.functions.hh>

// C++ Headers
#include <cctype>
#include <cstddef>
#include <cstring>
#include <string>


namespace ObjexxFCL {


// Constants
char const SPACE( ' ' );


// Predicate


/// @brief char == char Case-Optionally?
bool
equal( char const c, char const d, bool const exact_case )
{
	if ( exact_case ) {
		return ( c == d );
	} else {
		return ( std::tolower( c ) == std::tolower( d ) );
	}
}


/// @brief char == char Case-Insensitively
bool
equali( char const c, char const d )
{
	return ( std::tolower( c ) == std::tolower( d ) );
}


/// @brief Character is Blank?
bool
is_blank( char const c )
{
	return ( c == SPACE );
}


/// @brief Character is Not Blank?
bool
not_blank( char const c )
{
	return ( c != SPACE );
}


/// @brief Character is in a string?
bool
is_any_of( char const c, std::string const & s )
{
	return ( s.find( c ) != std::string::npos );
}


/// @brief Character is in a cstring?
bool
is_any_of( char const c, c_cstring const s )
{
	for ( std::size_t i = 0, e = std::strlen( s ); i < e; ++i ) {
		if ( c == s[ i ] ) return true;
	}
	return false; // No matches
}


// Modifier


/// @brief Lowercase a Character
char &
lowercase( char & c )
{
	c = std::tolower( c );
	return c;
}


/// @brief Uppercase a Character
char &
uppercase( char & c )
{
	c = std::toupper( c );
	return c;
}


// Generator


/// @brief Lowercased Copy of a Character
char
lowercased( char const c )
{
	return std::tolower( c );
}


/// @brief Uppercased Copy of a Character
char
uppercased( char const c )
{
	return std::toupper( c );
}


} // namespace ObjexxFCL
