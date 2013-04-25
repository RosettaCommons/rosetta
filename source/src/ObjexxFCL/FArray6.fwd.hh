#ifndef INCLUDED_ObjexxFCL_FArray6_fwd_hh
#define INCLUDED_ObjexxFCL_FArray6_fwd_hh


// FArray6 Forward Declarations
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
#include <cstddef>
#include <string>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray6;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray6< bool >                FArray6_bool;
typedef  FArray6< byte >                FArray6_byte;
typedef  FArray6< sbyte >               FArray6_sbyte;
typedef  FArray6< ubyte >               FArray6_ubyte;
typedef  FArray6< short int >           FArray6_short;
typedef  FArray6< int >                 FArray6_int;
typedef  FArray6< long int >            FArray6_long;
typedef  FArray6< unsigned short int >  FArray6_ushort;
typedef  FArray6< unsigned int >        FArray6_uint;
typedef  FArray6< unsigned long int >   FArray6_ulong;
typedef  FArray6< std::size_t >         FArray6_size_t;
typedef  FArray6< std::size_t >         FArray6_size;
typedef  FArray6< float >               FArray6_float;
typedef  FArray6< double >              FArray6_double;
typedef  FArray6< long double >         FArray6_longdouble;
typedef  FArray6< char >                FArray6_char;
typedef  FArray6< unsigned char >       FArray6_uchar;
typedef  FArray6< signed char >         FArray6_schar;
typedef  FArray6< std::string >         FArray6_string;
typedef  FArray6< Fstring >             FArray6_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray6_fwd_HH
