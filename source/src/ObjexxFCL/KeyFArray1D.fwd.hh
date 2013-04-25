#ifndef INCLUDED_ObjexxFCL_KeyFArray1D_fwd_hh
#define INCLUDED_ObjexxFCL_KeyFArray1D_fwd_hh


// KeyFArray1D Forward Declarations
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
template< typename > class KeyFArray1D;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  KeyFArray1D< bool >                KeyFArray1D_bool;
typedef  KeyFArray1D< byte >                KeyFArray1D_byte;
typedef  KeyFArray1D< sbyte >               KeyFArray1D_sbyte;
typedef  KeyFArray1D< ubyte >               KeyFArray1D_ubyte;
typedef  KeyFArray1D< short int >           KeyFArray1D_short;
typedef  KeyFArray1D< int >                 KeyFArray1D_int;
typedef  KeyFArray1D< long int >            KeyFArray1D_long;
typedef  KeyFArray1D< unsigned short int >  KeyFArray1D_ushort;
typedef  KeyFArray1D< unsigned int >        KeyFArray1D_uint;
typedef  KeyFArray1D< unsigned long int >   KeyFArray1D_ulong;
typedef  KeyFArray1D< std::size_t >         KeyFArray1D_size_t;
typedef  KeyFArray1D< std::size_t >         KeyFArray1D_size;
typedef  KeyFArray1D< float >               KeyFArray1D_float;
typedef  KeyFArray1D< double >              KeyFArray1D_double;
typedef  KeyFArray1D< long double >         KeyFArray1D_longdouble;
typedef  KeyFArray1D< char >                KeyFArray1D_char;
typedef  KeyFArray1D< unsigned char >       KeyFArray1D_uchar;
typedef  KeyFArray1D< signed char >         KeyFArray1D_schar;
typedef  KeyFArray1D< std::string >         KeyFArray1D_string;
typedef  KeyFArray1D< Fstring >             KeyFArray1D_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_KeyFArray1D_fwd_HH
