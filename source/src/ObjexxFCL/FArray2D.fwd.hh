#ifndef INCLUDED_ObjexxFCL_FArray2D_fwd_hh
#define INCLUDED_ObjexxFCL_FArray2D_fwd_hh


// FArray2D Forward Declarations
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
template< typename > class FArray2D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray2D< bool >                FArray2D_bool;
typedef  FArray2D< byte >                FArray2D_byte;
typedef  FArray2D< sbyte >               FArray2D_sbyte;
typedef  FArray2D< ubyte >               FArray2D_ubyte;
typedef  FArray2D< short int >           FArray2D_short;
typedef  FArray2D< int >                 FArray2D_int;
typedef  FArray2D< long int >            FArray2D_long;
typedef  FArray2D< unsigned short int >  FArray2D_ushort;
typedef  FArray2D< unsigned int >        FArray2D_uint;
typedef  FArray2D< unsigned long int >   FArray2D_ulong;
typedef  FArray2D< std::size_t >         FArray2D_size_t;
typedef  FArray2D< std::size_t >         FArray2D_size;
typedef  FArray2D< float >               FArray2D_float;
typedef  FArray2D< double >              FArray2D_double;
typedef  FArray2D< long double >         FArray2D_longdouble;
typedef  FArray2D< char >                FArray2D_char;
typedef  FArray2D< unsigned char >       FArray2D_uchar;
typedef  FArray2D< signed char >         FArray2D_schar;
typedef  FArray2D< std::string >         FArray2D_string;
typedef  FArray2D< Fstring >             FArray2D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray2D_fwd_HH
