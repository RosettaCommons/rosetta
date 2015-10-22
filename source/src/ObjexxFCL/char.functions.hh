#ifndef INCLUDED_ObjexxFCL_char_functions_hh
#define INCLUDED_ObjexxFCL_char_functions_hh


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


// C++ Headers
#ifdef WIN_PYROSETTA
#include <string>
#else
#include <iosfwd>
#endif

namespace ObjexxFCL {


// Types
typedef  char const *  c_cstring;


// Predicate


/// @brief char == char Case-Optionally?
bool
equal( char const c, char const d, bool const exact_case = true );


/// @brief char == char Case-Insensitively
bool
equali( char const c, char const d );


/// @brief Character is Blank?
bool
is_blank( char const c );


/// @brief Character is Not Blank?
bool
not_blank( char const c );


/// @brief Character is in a string?
bool
is_any_of( char const c, std::string const & s );


/// @brief Character is in a cstring?
bool
is_any_of( char const c, c_cstring const s );


/// @brief ASCII Lexical >= Comparison
inline
bool
lge( char const s, char const t )
{
	return ( s >= t );
}


/// @brief ASCII Lexical < Comparison
inline
bool
lgt( char const s, char const t )
{
	return ( s > t );
}


/// @brief ASCII Lexical <= Comparison
inline
bool
lle( char const s, char const t )
{
	return ( s <= t );
}


/// @brief ASCII Lexical < Comparison
inline
bool
llt( char const s, char const t )
{
	return ( s < t );
}


// Integer Conversion


/// @brief Integer Value of a Given One-Character Fstring
inline
int
ICHAR( char const s )
{
	return static_cast< int >( s );
}


/// @brief ASCII Integer Value for a Given One-Character Fstring
inline
int
IACHAR( char const s )
{
	return static_cast< int >( s );
}


// Modifier


/// @brief Lowercase a Character
char &
lowercase( char & c );


/// @brief Uppercase a Character
char &
uppercase( char & c );


// Generator


/// @brief Space Character
inline
char
space()
{
	return ' ';
}


/// @brief Lowercased Copy of a Character
char
lowercased( char const c );


/// @brief Uppercased Copy of a Character
char
uppercased( char const c );


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_char_functions_HH
