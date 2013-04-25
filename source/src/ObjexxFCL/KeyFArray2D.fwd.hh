#ifndef INCLUDED_ObjexxFCL_KeyFArray2D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray2D_fwd_hh


// KeyFArray2D Forward Declarations
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
template< typename > class KeyFArray2D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray2D< bool >                KeyFArray2D_bool;
typedef  KeyFArray2D< byte >                KeyFArray2D_byte;
typedef  KeyFArray2D< sbyte >               KeyFArray2D_sbyte;
typedef  KeyFArray2D< ubyte >               KeyFArray2D_ubyte;
typedef  KeyFArray2D< short int >           KeyFArray2D_short;
typedef  KeyFArray2D< int >                 KeyFArray2D_int;
typedef  KeyFArray2D< long int >            KeyFArray2D_long;
typedef  KeyFArray2D< unsigned short int >  KeyFArray2D_ushort;
typedef  KeyFArray2D< unsigned int >        KeyFArray2D_uint;
typedef  KeyFArray2D< unsigned long int >   KeyFArray2D_ulong;
typedef  KeyFArray2D< std::size_t >         KeyFArray2D_size_t;
typedef  KeyFArray2D< std::size_t >         KeyFArray2D_size;
typedef  KeyFArray2D< float >               KeyFArray2D_float;
typedef  KeyFArray2D< double >              KeyFArray2D_double;
typedef  KeyFArray2D< long double >         KeyFArray2D_longdouble;
typedef  KeyFArray2D< char >                KeyFArray2D_char;
typedef  KeyFArray2D< unsigned char >       KeyFArray2D_uchar;
typedef  KeyFArray2D< signed char >         KeyFArray2D_schar;
typedef  KeyFArray2D< std::string >         KeyFArray2D_string;
typedef  KeyFArray2D< Fstring >             KeyFArray2D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray2D_fwd_HH
