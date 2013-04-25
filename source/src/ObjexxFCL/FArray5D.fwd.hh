#ifndef INCLUDED_ObjexxFCL_FArray5D_fwd_hh
#define INCLUDED_ObjexxFCL_FArray5D_fwd_hh


// FArray5D Forward Declarations
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
template< typename > class FArray5D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray5D< bool >                FArray5D_bool;
typedef  FArray5D< byte >                FArray5D_byte;
typedef  FArray5D< sbyte >               FArray5D_sbyte;
typedef  FArray5D< ubyte >               FArray5D_ubyte;
typedef  FArray5D< short int >           FArray5D_short;
typedef  FArray5D< int >                 FArray5D_int;
typedef  FArray5D< long int >            FArray5D_long;
typedef  FArray5D< unsigned short int >  FArray5D_ushort;
typedef  FArray5D< unsigned int >        FArray5D_uint;
typedef  FArray5D< unsigned long int >   FArray5D_ulong;
typedef  FArray5D< std::size_t >         FArray5D_size_t;
typedef  FArray5D< std::size_t >         FArray5D_size;
typedef  FArray5D< float >               FArray5D_float;
typedef  FArray5D< double >              FArray5D_double;
typedef  FArray5D< long double >         FArray5D_longdouble;
typedef  FArray5D< char >                FArray5D_char;
typedef  FArray5D< unsigned char >       FArray5D_uchar;
typedef  FArray5D< signed char >         FArray5D_schar;
typedef  FArray5D< std::string >         FArray5D_string;
typedef  FArray5D< Fstring >             FArray5D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray5D_fwd_HH
