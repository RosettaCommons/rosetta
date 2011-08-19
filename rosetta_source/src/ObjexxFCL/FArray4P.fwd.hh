#ifndef INCLUDED_ObjexxFCL_FArray4P_fwd_hh
#define INCLUDED_ObjexxFCL_FArray4P_fwd_hh


// FArray4P Forward Declarations
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
template< typename > class FArray4P;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray4P< bool >                FArray4P_bool;
typedef  FArray4P< byte >                FArray4P_byte;
typedef  FArray4P< sbyte >               FArray4P_sbyte;
typedef  FArray4P< ubyte >               FArray4P_ubyte;
typedef  FArray4P< short int >           FArray4P_short;
typedef  FArray4P< int >                 FArray4P_int;
typedef  FArray4P< long int >            FArray4P_long;
typedef  FArray4P< unsigned short int >  FArray4P_ushort;
typedef  FArray4P< unsigned int >        FArray4P_uint;
typedef  FArray4P< unsigned long int >   FArray4P_ulong;
typedef  FArray4P< std::size_t >         FArray4P_size_t;
typedef  FArray4P< std::size_t >         FArray4P_size;
typedef  FArray4P< float >               FArray4P_float;
typedef  FArray4P< double >              FArray4P_double;
typedef  FArray4P< long double >         FArray4P_longdouble;
typedef  FArray4P< char >                FArray4P_char;
typedef  FArray4P< unsigned char >       FArray4P_uchar;
typedef  FArray4P< signed char >         FArray4P_schar;
typedef  FArray4P< std::string >         FArray4P_string;
typedef  FArray4P< Fstring >             FArray4P_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray4P_fwd_HH
