#ifndef INCLUDED_ObjexxFCL_FArray1_fwd_hh
#define INCLUDED_ObjexxFCL_FArray1_fwd_hh


// FArray1 Forward Declarations
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
template< typename > class FArray1;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray1< bool >                FArray1_bool;
typedef  FArray1< byte >                FArray1_byte;
typedef  FArray1< sbyte >               FArray1_sbyte;
typedef  FArray1< ubyte >               FArray1_ubyte;
typedef  FArray1< short int >           FArray1_short;
typedef  FArray1< int >                 FArray1_int;
typedef  FArray1< long int >            FArray1_long;
typedef  FArray1< unsigned short int >  FArray1_ushort;
typedef  FArray1< unsigned int >        FArray1_uint;
typedef  FArray1< unsigned long int >   FArray1_ulong;
typedef  FArray1< std::size_t >         FArray1_size_t;
typedef  FArray1< std::size_t >         FArray1_size;
typedef  FArray1< float >               FArray1_float;
typedef  FArray1< double >              FArray1_double;
typedef  FArray1< long double >         FArray1_longdouble;
typedef  FArray1< char >                FArray1_char;
typedef  FArray1< unsigned char >       FArray1_uchar;
typedef  FArray1< signed char >         FArray1_schar;
typedef  FArray1< std::string >         FArray1_string;
typedef  FArray1< Fstring >             FArray1_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1_fwd_HH
