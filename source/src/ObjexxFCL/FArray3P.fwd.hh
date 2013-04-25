#ifndef INCLUDED_ObjexxFCL_FArray3P_fwd_hh
#define INCLUDED_ObjexxFCL_FArray3P_fwd_hh


// FArray3P Forward Declarations
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
template< typename > class FArray3P;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray3P< bool >                FArray3P_bool;
typedef  FArray3P< byte >                FArray3P_byte;
typedef  FArray3P< sbyte >               FArray3P_sbyte;
typedef  FArray3P< ubyte >               FArray3P_ubyte;
typedef  FArray3P< short int >           FArray3P_short;
typedef  FArray3P< int >                 FArray3P_int;
typedef  FArray3P< long int >            FArray3P_long;
typedef  FArray3P< unsigned short int >  FArray3P_ushort;
typedef  FArray3P< unsigned int >        FArray3P_uint;
typedef  FArray3P< unsigned long int >   FArray3P_ulong;
typedef  FArray3P< std::size_t >         FArray3P_size_t;
typedef  FArray3P< std::size_t >         FArray3P_size;
typedef  FArray3P< float >               FArray3P_float;
typedef  FArray3P< double >              FArray3P_double;
typedef  FArray3P< long double >         FArray3P_longdouble;
typedef  FArray3P< char >                FArray3P_char;
typedef  FArray3P< unsigned char >       FArray3P_uchar;
typedef  FArray3P< signed char >         FArray3P_schar;
typedef  FArray3P< std::string >         FArray3P_string;
typedef  FArray3P< Fstring >             FArray3P_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray3P_fwd_HH
