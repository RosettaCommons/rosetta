#ifndef INCLUDED_ObjexxFCL_CArray_fwd_hh
#define INCLUDED_ObjexxFCL_CArray_fwd_hh


// CArray Forward Declarations
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
template< typename > class CArray;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  CArray< bool >                CArray_bool;
typedef  CArray< byte >                CArray_byte;
typedef  CArray< sbyte >               CArray_sbyte;
typedef  CArray< ubyte >               CArray_ubyte;
typedef  CArray< short int >           CArray_short;
typedef  CArray< int >                 CArray_int;
typedef  CArray< long int >            CArray_long;
typedef  CArray< unsigned short int >  CArray_ushort;
typedef  CArray< unsigned int >        CArray_uint;
typedef  CArray< unsigned long int >   CArray_ulong;
typedef  CArray< std::size_t >         CArray_size_t;
typedef  CArray< std::size_t >         CArray_size;
typedef  CArray< float >               CArray_float;
typedef  CArray< double >              CArray_double;
typedef  CArray< long double >         CArray_longdouble;
typedef  CArray< char >                CArray_char;
typedef  CArray< unsigned char >       CArray_uchar;
typedef  CArray< signed char >         CArray_schar;
typedef  CArray< std::string >         CArray_string;
typedef  CArray< Fstring >             CArray_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_CArray_fwd_HH
