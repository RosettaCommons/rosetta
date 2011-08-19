#ifndef INCLUDED_ObjexxFCL_KeyFArray3D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray3D_fwd_hh


// KeyFArray3D Forward Declarations
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
template< typename > class KeyFArray3D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray3D< bool >                KeyFArray3D_bool;
typedef  KeyFArray3D< byte >                KeyFArray3D_byte;
typedef  KeyFArray3D< sbyte >               KeyFArray3D_sbyte;
typedef  KeyFArray3D< ubyte >               KeyFArray3D_ubyte;
typedef  KeyFArray3D< short int >           KeyFArray3D_short;
typedef  KeyFArray3D< int >                 KeyFArray3D_int;
typedef  KeyFArray3D< long int >            KeyFArray3D_long;
typedef  KeyFArray3D< unsigned short int >  KeyFArray3D_ushort;
typedef  KeyFArray3D< unsigned int >        KeyFArray3D_uint;
typedef  KeyFArray3D< unsigned long int >   KeyFArray3D_ulong;
typedef  KeyFArray3D< std::size_t >         KeyFArray3D_size_t;
typedef  KeyFArray3D< std::size_t >         KeyFArray3D_size;
typedef  KeyFArray3D< float >               KeyFArray3D_float;
typedef  KeyFArray3D< double >              KeyFArray3D_double;
typedef  KeyFArray3D< long double >         KeyFArray3D_longdouble;
typedef  KeyFArray3D< char >                KeyFArray3D_char;
typedef  KeyFArray3D< unsigned char >       KeyFArray3D_uchar;
typedef  KeyFArray3D< signed char >         KeyFArray3D_schar;
typedef  KeyFArray3D< std::string >         KeyFArray3D_string;
typedef  KeyFArray3D< Fstring >             KeyFArray3D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray3D_fwd_HH
