#ifndef INCLUDED_ObjexxFCL_KeyFArray6D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray6D_fwd_hh


// KeyFArray6D Forward Declarations
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
template< typename > class KeyFArray6D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray6D< bool >                KeyFArray6D_bool;
typedef  KeyFArray6D< byte >                KeyFArray6D_byte;
typedef  KeyFArray6D< sbyte >               KeyFArray6D_sbyte;
typedef  KeyFArray6D< ubyte >               KeyFArray6D_ubyte;
typedef  KeyFArray6D< short int >           KeyFArray6D_short;
typedef  KeyFArray6D< int >                 KeyFArray6D_int;
typedef  KeyFArray6D< long int >            KeyFArray6D_long;
typedef  KeyFArray6D< unsigned short int >  KeyFArray6D_ushort;
typedef  KeyFArray6D< unsigned int >        KeyFArray6D_uint;
typedef  KeyFArray6D< unsigned long int >   KeyFArray6D_ulong;
typedef  KeyFArray6D< std::size_t >         KeyFArray6D_size_t;
typedef  KeyFArray6D< std::size_t >         KeyFArray6D_size;
typedef  KeyFArray6D< float >               KeyFArray6D_float;
typedef  KeyFArray6D< double >              KeyFArray6D_double;
typedef  KeyFArray6D< long double >         KeyFArray6D_longdouble;
typedef  KeyFArray6D< char >                KeyFArray6D_char;
typedef  KeyFArray6D< unsigned char >       KeyFArray6D_uchar;
typedef  KeyFArray6D< signed char >         KeyFArray6D_schar;
typedef  KeyFArray6D< std::string >         KeyFArray6D_string;
typedef  KeyFArray6D< Fstring >             KeyFArray6D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray6D_fwd_HH
