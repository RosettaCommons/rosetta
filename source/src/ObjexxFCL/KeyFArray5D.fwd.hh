#ifndef INCLUDED_ObjexxFCL_KeyFArray5D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray5D_fwd_hh


// KeyFArray5D Forward Declarations
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
template< typename > class KeyFArray5D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray5D< bool >                KeyFArray5D_bool;
typedef  KeyFArray5D< byte >                KeyFArray5D_byte;
typedef  KeyFArray5D< sbyte >               KeyFArray5D_sbyte;
typedef  KeyFArray5D< ubyte >               KeyFArray5D_ubyte;
typedef  KeyFArray5D< short int >           KeyFArray5D_short;
typedef  KeyFArray5D< int >                 KeyFArray5D_int;
typedef  KeyFArray5D< long int >            KeyFArray5D_long;
typedef  KeyFArray5D< unsigned short int >  KeyFArray5D_ushort;
typedef  KeyFArray5D< unsigned int >        KeyFArray5D_uint;
typedef  KeyFArray5D< unsigned long int >   KeyFArray5D_ulong;
typedef  KeyFArray5D< std::size_t >         KeyFArray5D_size_t;
typedef  KeyFArray5D< std::size_t >         KeyFArray5D_size;
typedef  KeyFArray5D< float >               KeyFArray5D_float;
typedef  KeyFArray5D< double >              KeyFArray5D_double;
typedef  KeyFArray5D< long double >         KeyFArray5D_longdouble;
typedef  KeyFArray5D< char >                KeyFArray5D_char;
typedef  KeyFArray5D< unsigned char >       KeyFArray5D_uchar;
typedef  KeyFArray5D< signed char >         KeyFArray5D_schar;
typedef  KeyFArray5D< std::string >         KeyFArray5D_string;
typedef  KeyFArray5D< Fstring >             KeyFArray5D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray5D_fwd_HH
