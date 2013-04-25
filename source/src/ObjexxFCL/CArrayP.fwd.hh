#ifndef INCLUDED_ObjexxFCL_CArrayP_fwd_hh
#define INCLUDED_ObjexxFCL_CArrayP_fwd_hh


// CArrayP Forward Declarations
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
template< typename > class CArrayP;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  CArrayP< bool >                CArrayP_bool;
typedef  CArrayP< byte >                CArrayP_byte;
typedef  CArrayP< sbyte >               CArrayP_sbyte;
typedef  CArrayP< ubyte >               CArrayP_ubyte;
typedef  CArrayP< short int >           CArrayP_short;
typedef  CArrayP< int >                 CArrayP_int;
typedef  CArrayP< long int >            CArrayP_long;
typedef  CArrayP< unsigned short int >  CArrayP_ushort;
typedef  CArrayP< unsigned int >        CArrayP_uint;
typedef  CArrayP< unsigned long int >   CArrayP_ulong;
typedef  CArrayP< std::size_t >         CArrayP_size_t;
typedef  CArrayP< std::size_t >         CArrayP_size;
typedef  CArrayP< float >               CArrayP_float;
typedef  CArrayP< double >              CArrayP_double;
typedef  CArrayP< long double >         CArrayP_longdouble;
typedef  CArrayP< char >                CArrayP_char;
typedef  CArrayP< unsigned char >       CArrayP_uchar;
typedef  CArrayP< signed char >         CArrayP_schar;
typedef  CArrayP< std::string >         CArrayP_string;
typedef  CArrayP< Fstring >             CArrayP_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_CArrayP_fwd_HH
