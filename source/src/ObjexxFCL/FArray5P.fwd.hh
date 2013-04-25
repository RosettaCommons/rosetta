#ifndef INCLUDED_ObjexxFCL_FArray5P_fwd_hh
#define INCLUDED_ObjexxFCL_FArray5P_fwd_hh


// FArray5P Forward Declarations
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
template< typename > class FArray5P;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray5P< bool >                FArray5P_bool;
typedef  FArray5P< byte >                FArray5P_byte;
typedef  FArray5P< sbyte >               FArray5P_sbyte;
typedef  FArray5P< ubyte >               FArray5P_ubyte;
typedef  FArray5P< short int >           FArray5P_short;
typedef  FArray5P< int >                 FArray5P_int;
typedef  FArray5P< long int >            FArray5P_long;
typedef  FArray5P< unsigned short int >  FArray5P_ushort;
typedef  FArray5P< unsigned int >        FArray5P_uint;
typedef  FArray5P< unsigned long int >   FArray5P_ulong;
typedef  FArray5P< std::size_t >         FArray5P_size_t;
typedef  FArray5P< std::size_t >         FArray5P_size;
typedef  FArray5P< float >               FArray5P_float;
typedef  FArray5P< double >              FArray5P_double;
typedef  FArray5P< long double >         FArray5P_longdouble;
typedef  FArray5P< char >                FArray5P_char;
typedef  FArray5P< unsigned char >       FArray5P_uchar;
typedef  FArray5P< signed char >         FArray5P_schar;
typedef  FArray5P< std::string >         FArray5P_string;
typedef  FArray5P< Fstring >             FArray5P_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray5P_fwd_HH
