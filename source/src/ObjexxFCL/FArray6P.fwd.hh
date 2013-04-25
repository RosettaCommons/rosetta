#ifndef INCLUDED_ObjexxFCL_FArray6P_fwd_hh
#define INCLUDED_ObjexxFCL_FArray6P_fwd_hh


// FArray6P Forward Declarations
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
template< typename > class FArray6P;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray6P< bool >                FArray6P_bool;
typedef  FArray6P< byte >                FArray6P_byte;
typedef  FArray6P< sbyte >               FArray6P_sbyte;
typedef  FArray6P< ubyte >               FArray6P_ubyte;
typedef  FArray6P< short int >           FArray6P_short;
typedef  FArray6P< int >                 FArray6P_int;
typedef  FArray6P< long int >            FArray6P_long;
typedef  FArray6P< unsigned short int >  FArray6P_ushort;
typedef  FArray6P< unsigned int >        FArray6P_uint;
typedef  FArray6P< unsigned long int >   FArray6P_ulong;
typedef  FArray6P< std::size_t >         FArray6P_size_t;
typedef  FArray6P< std::size_t >         FArray6P_size;
typedef  FArray6P< float >               FArray6P_float;
typedef  FArray6P< double >              FArray6P_double;
typedef  FArray6P< long double >         FArray6P_longdouble;
typedef  FArray6P< char >                FArray6P_char;
typedef  FArray6P< unsigned char >       FArray6P_uchar;
typedef  FArray6P< signed char >         FArray6P_schar;
typedef  FArray6P< std::string >         FArray6P_string;
typedef  FArray6P< Fstring >             FArray6P_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray6P_fwd_HH
