#ifndef INCLUDED_ObjexxFCL_FArray_fwd_hh
#define INCLUDED_ObjexxFCL_FArray_fwd_hh


// FArray Forward Declarations
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
template< typename > class FArray;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray< bool >                FArray_bool;
typedef  FArray< byte >                FArray_byte;
typedef  FArray< sbyte >               FArray_sbyte;
typedef  FArray< ubyte >               FArray_ubyte;
typedef  FArray< short int >           FArray_short;
typedef  FArray< int >                 FArray_int;
typedef  FArray< long int >            FArray_long;
typedef  FArray< unsigned short int >  FArray_ushort;
typedef  FArray< unsigned int >        FArray_uint;
typedef  FArray< unsigned long int >   FArray_ulong;
typedef  FArray< std::size_t >         FArray_size_t;
typedef  FArray< std::size_t >         FArray_size;
typedef  FArray< float >               FArray_float;
typedef  FArray< double >              FArray_double;
typedef  FArray< long double >         FArray_longdouble;
typedef  FArray< char >                FArray_char;
typedef  FArray< unsigned char >       FArray_uchar;
typedef  FArray< signed char >         FArray_schar;
typedef  FArray< std::string >         FArray_string;
typedef  FArray< Fstring >             FArray_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray_fwd_HH
