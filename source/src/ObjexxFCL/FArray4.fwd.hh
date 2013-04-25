#ifndef INCLUDED_ObjexxFCL_FArray4_fwd_hh
#define INCLUDED_ObjexxFCL_FArray4_fwd_hh


// FArray4 Forward Declarations
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
template< typename > class FArray4;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray4< bool >                FArray4_bool;
typedef  FArray4< byte >                FArray4_byte;
typedef  FArray4< sbyte >               FArray4_sbyte;
typedef  FArray4< ubyte >               FArray4_ubyte;
typedef  FArray4< short int >           FArray4_short;
typedef  FArray4< int >                 FArray4_int;
typedef  FArray4< long int >            FArray4_long;
typedef  FArray4< unsigned short int >  FArray4_ushort;
typedef  FArray4< unsigned int >        FArray4_uint;
typedef  FArray4< unsigned long int >   FArray4_ulong;
typedef  FArray4< std::size_t >         FArray4_size_t;
typedef  FArray4< std::size_t >         FArray4_size;
typedef  FArray4< float >               FArray4_float;
typedef  FArray4< double >              FArray4_double;
typedef  FArray4< long double >         FArray4_longdouble;
typedef  FArray4< char >                FArray4_char;
typedef  FArray4< unsigned char >       FArray4_uchar;
typedef  FArray4< signed char >         FArray4_schar;
typedef  FArray4< std::string >         FArray4_string;
typedef  FArray4< Fstring >             FArray4_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray4_fwd_HH
