#ifndef INCLUDED_ObjexxFCL_FArray1D_fwd_hh
#define INCLUDED_ObjexxFCL_FArray1D_fwd_hh


// FArray1D Forward Declarations
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
template< typename > class FArray1D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray1D< bool >                FArray1D_bool;
typedef  FArray1D< byte >                FArray1D_byte;
typedef  FArray1D< sbyte >               FArray1D_sbyte;
typedef  FArray1D< ubyte >               FArray1D_ubyte;
typedef  FArray1D< short int >           FArray1D_short;
typedef  FArray1D< int >                 FArray1D_int;
typedef  FArray1D< long int >            FArray1D_long;
typedef  FArray1D< unsigned short int >  FArray1D_ushort;
typedef  FArray1D< unsigned int >        FArray1D_uint;
typedef  FArray1D< unsigned long int >   FArray1D_ulong;
typedef  FArray1D< std::size_t >         FArray1D_size_t;
typedef  FArray1D< std::size_t >         FArray1D_size;
typedef  FArray1D< float >               FArray1D_float;
typedef  FArray1D< double >              FArray1D_double;
typedef  FArray1D< long double >         FArray1D_longdouble;
typedef  FArray1D< char >                FArray1D_char;
typedef  FArray1D< unsigned char >       FArray1D_uchar;
typedef  FArray1D< signed char >         FArray1D_schar;
typedef  FArray1D< std::string >         FArray1D_string;
typedef  FArray1D< Fstring >             FArray1D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1D_fwd_HH
