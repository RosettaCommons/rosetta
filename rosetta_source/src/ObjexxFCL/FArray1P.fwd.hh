#ifndef INCLUDED_ObjexxFCL_FArray1P_fwd_hh
#define INCLUDED_ObjexxFCL_FArray1P_fwd_hh


// FArray1P Forward Declarations
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
template< typename > class FArray1P;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray1P< bool >                FArray1P_bool;
typedef  FArray1P< byte >                FArray1P_byte;
typedef  FArray1P< sbyte >               FArray1P_sbyte;
typedef  FArray1P< ubyte >               FArray1P_ubyte;
typedef  FArray1P< short int >           FArray1P_short;
typedef  FArray1P< int >                 FArray1P_int;
typedef  FArray1P< long int >            FArray1P_long;
typedef  FArray1P< unsigned short int >  FArray1P_ushort;
typedef  FArray1P< unsigned int >        FArray1P_uint;
typedef  FArray1P< unsigned long int >   FArray1P_ulong;
typedef  FArray1P< std::size_t >         FArray1P_size_t;
typedef  FArray1P< std::size_t >         FArray1P_size;
typedef  FArray1P< float >               FArray1P_float;
typedef  FArray1P< double >              FArray1P_double;
typedef  FArray1P< long double >         FArray1P_longdouble;
typedef  FArray1P< char >                FArray1P_char;
typedef  FArray1P< unsigned char >       FArray1P_uchar;
typedef  FArray1P< signed char >         FArray1P_schar;
typedef  FArray1P< std::string >         FArray1P_string;
typedef  FArray1P< Fstring >             FArray1P_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1P_fwd_HH
