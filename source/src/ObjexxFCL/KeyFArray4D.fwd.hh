#ifndef INCLUDED_ObjexxFCL_KeyFArray4D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray4D_fwd_hh


// KeyFArray4D Forward Declarations
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
template< typename > class KeyFArray4D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray4D< bool >                KeyFArray4D_bool;
typedef  KeyFArray4D< byte >                KeyFArray4D_byte;
typedef  KeyFArray4D< sbyte >               KeyFArray4D_sbyte;
typedef  KeyFArray4D< ubyte >               KeyFArray4D_ubyte;
typedef  KeyFArray4D< short int >           KeyFArray4D_short;
typedef  KeyFArray4D< int >                 KeyFArray4D_int;
typedef  KeyFArray4D< long int >            KeyFArray4D_long;
typedef  KeyFArray4D< unsigned short int >  KeyFArray4D_ushort;
typedef  KeyFArray4D< unsigned int >        KeyFArray4D_uint;
typedef  KeyFArray4D< unsigned long int >   KeyFArray4D_ulong;
typedef  KeyFArray4D< std::size_t >         KeyFArray4D_size_t;
typedef  KeyFArray4D< std::size_t >         KeyFArray4D_size;
typedef  KeyFArray4D< float >               KeyFArray4D_float;
typedef  KeyFArray4D< double >              KeyFArray4D_double;
typedef  KeyFArray4D< long double >         KeyFArray4D_longdouble;
typedef  KeyFArray4D< char >                KeyFArray4D_char;
typedef  KeyFArray4D< unsigned char >       KeyFArray4D_uchar;
typedef  KeyFArray4D< signed char >         KeyFArray4D_schar;
typedef  KeyFArray4D< std::string >         KeyFArray4D_string;
typedef  KeyFArray4D< Fstring >             KeyFArray4D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray4D_fwd_HH
