#ifndef INCLUDED_ObjexxFCL_FArray2_fwd_hh
#define INCLUDED_ObjexxFCL_FArray2_fwd_hh


// FArray2 Forward Declarations
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
template< typename > class FArray2;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray2< bool >                FArray2_bool;
typedef  FArray2< byte >                FArray2_byte;
typedef  FArray2< sbyte >               FArray2_sbyte;
typedef  FArray2< ubyte >               FArray2_ubyte;
typedef  FArray2< short int >           FArray2_short;
typedef  FArray2< int >                 FArray2_int;
typedef  FArray2< long int >            FArray2_long;
typedef  FArray2< unsigned short int >  FArray2_ushort;
typedef  FArray2< unsigned int >        FArray2_uint;
typedef  FArray2< unsigned long int >   FArray2_ulong;
typedef  FArray2< std::size_t >         FArray2_size_t;
typedef  FArray2< std::size_t >         FArray2_size;
typedef  FArray2< float >               FArray2_float;
typedef  FArray2< double >              FArray2_double;
typedef  FArray2< long double >         FArray2_longdouble;
typedef  FArray2< char >                FArray2_char;
typedef  FArray2< unsigned char >       FArray2_uchar;
typedef  FArray2< signed char >         FArray2_schar;
typedef  FArray2< std::string >         FArray2_string;
typedef  FArray2< Fstring >             FArray2_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray2_fwd_HH
