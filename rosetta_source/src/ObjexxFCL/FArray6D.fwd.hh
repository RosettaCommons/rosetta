#ifndef INCLUDED_ObjexxFCL_FArray6D_fwd_hh
#define INCLUDED_ObjexxFCL_FArray6D_fwd_hh


// FArray6D Forward Declarations
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
template< typename > class FArray6D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray6D< bool >                FArray6D_bool;
typedef  FArray6D< byte >                FArray6D_byte;
typedef  FArray6D< sbyte >               FArray6D_sbyte;
typedef  FArray6D< ubyte >               FArray6D_ubyte;
typedef  FArray6D< short int >           FArray6D_short;
typedef  FArray6D< int >                 FArray6D_int;
typedef  FArray6D< long int >            FArray6D_long;
typedef  FArray6D< unsigned short int >  FArray6D_ushort;
typedef  FArray6D< unsigned int >        FArray6D_uint;
typedef  FArray6D< unsigned long int >   FArray6D_ulong;
typedef  FArray6D< std::size_t >         FArray6D_size_t;
typedef  FArray6D< std::size_t >         FArray6D_size;
typedef  FArray6D< float >               FArray6D_float;
typedef  FArray6D< double >              FArray6D_double;
typedef  FArray6D< long double >         FArray6D_longdouble;
typedef  FArray6D< char >                FArray6D_char;
typedef  FArray6D< unsigned char >       FArray6D_uchar;
typedef  FArray6D< signed char >         FArray6D_schar;
typedef  FArray6D< std::string >         FArray6D_string;
typedef  FArray6D< Fstring >             FArray6D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray6D_fwd_HH
