#ifndef INCLUDED_ObjexxFCL_FArray1A_fwd_hh
#define INCLUDED_ObjexxFCL_FArray1A_fwd_hh


// FArray1A Forward Declarations
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
template< typename > class FArray1A;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray1A< bool >                FArray1A_bool;
typedef  FArray1A< byte >                FArray1A_byte;
typedef  FArray1A< sbyte >               FArray1A_sbyte;
typedef  FArray1A< ubyte >               FArray1A_ubyte;
typedef  FArray1A< short int >           FArray1A_short;
typedef  FArray1A< int >                 FArray1A_int;
typedef  FArray1A< long int >            FArray1A_long;
typedef  FArray1A< unsigned short int >  FArray1A_ushort;
typedef  FArray1A< unsigned int >        FArray1A_uint;
typedef  FArray1A< unsigned long int >   FArray1A_ulong;
typedef  FArray1A< std::size_t >         FArray1A_size_t;
typedef  FArray1A< std::size_t >         FArray1A_size;
typedef  FArray1A< float >               FArray1A_float;
typedef  FArray1A< double >              FArray1A_double;
typedef  FArray1A< long double >         FArray1A_longdouble;
typedef  FArray1A< char >                FArray1A_char;
typedef  FArray1A< unsigned char >       FArray1A_uchar;
typedef  FArray1A< signed char >         FArray1A_schar;
typedef  FArray1A< std::string >         FArray1A_string;
typedef  FArray1A< Fstring >             FArray1A_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1A_fwd_HH
