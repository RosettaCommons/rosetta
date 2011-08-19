#ifndef INCLUDED_ObjexxFCL_FArray3D_fwd_hh
#define INCLUDED_ObjexxFCL_FArray3D_fwd_hh


// FArray3D Forward Declarations
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
template< typename > class FArray3D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray3D< bool >                FArray3D_bool;
typedef  FArray3D< byte >                FArray3D_byte;
typedef  FArray3D< sbyte >               FArray3D_sbyte;
typedef  FArray3D< ubyte >               FArray3D_ubyte;
typedef  FArray3D< short int >           FArray3D_short;
typedef  FArray3D< int >                 FArray3D_int;
typedef  FArray3D< long int >            FArray3D_long;
typedef  FArray3D< unsigned short int >  FArray3D_ushort;
typedef  FArray3D< unsigned int >        FArray3D_uint;
typedef  FArray3D< unsigned long int >   FArray3D_ulong;
typedef  FArray3D< std::size_t >         FArray3D_size_t;
typedef  FArray3D< std::size_t >         FArray3D_size;
typedef  FArray3D< float >               FArray3D_float;
typedef  FArray3D< double >              FArray3D_double;
typedef  FArray3D< long double >         FArray3D_longdouble;
typedef  FArray3D< char >                FArray3D_char;
typedef  FArray3D< unsigned char >       FArray3D_uchar;
typedef  FArray3D< signed char >         FArray3D_schar;
typedef  FArray3D< std::string >         FArray3D_string;
typedef  FArray3D< Fstring >             FArray3D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray3D_fwd_HH
