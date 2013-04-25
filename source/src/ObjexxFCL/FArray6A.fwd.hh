#ifndef INCLUDED_ObjexxFCL_FArray6A_fwd_hh
#define INCLUDED_ObjexxFCL_FArray6A_fwd_hh


// FArray6A Forward Declarations
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
template< typename > class FArray6A;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray6A< bool >                FArray6A_bool;
typedef  FArray6A< byte >                FArray6A_byte;
typedef  FArray6A< sbyte >               FArray6A_sbyte;
typedef  FArray6A< ubyte >               FArray6A_ubyte;
typedef  FArray6A< short int >           FArray6A_short;
typedef  FArray6A< int >                 FArray6A_int;
typedef  FArray6A< long int >            FArray6A_long;
typedef  FArray6A< unsigned short int >  FArray6A_ushort;
typedef  FArray6A< unsigned int >        FArray6A_uint;
typedef  FArray6A< unsigned long int >   FArray6A_ulong;
typedef  FArray6A< std::size_t >         FArray6A_size_t;
typedef  FArray6A< std::size_t >         FArray6A_size;
typedef  FArray6A< float >               FArray6A_float;
typedef  FArray6A< double >              FArray6A_double;
typedef  FArray6A< long double >         FArray6A_longdouble;
typedef  FArray6A< char >                FArray6A_char;
typedef  FArray6A< unsigned char >       FArray6A_uchar;
typedef  FArray6A< signed char >         FArray6A_schar;
typedef  FArray6A< std::string >         FArray6A_string;
typedef  FArray6A< Fstring >             FArray6A_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray6A_fwd_HH
