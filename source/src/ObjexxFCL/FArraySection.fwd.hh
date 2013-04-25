#ifndef INCLUDED_ObjexxFCL_FArraySection_fwd_hh
#define INCLUDED_ObjexxFCL_FArraySection_fwd_hh


// FArraySection Forward Declarations
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
template< typename > class FArraySection;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArraySection< bool >                FArraySection_bool;
typedef  FArraySection< byte >                FArraySection_byte;
typedef  FArraySection< sbyte >               FArraySection_sbyte;
typedef  FArraySection< ubyte >               FArraySection_ubyte;
typedef  FArraySection< short int >           FArraySection_short;
typedef  FArraySection< int >                 FArraySection_int;
typedef  FArraySection< long int >            FArraySection_long;
typedef  FArraySection< unsigned short int >  FArraySection_ushort;
typedef  FArraySection< unsigned int >        FArraySection_uint;
typedef  FArraySection< unsigned long int >   FArraySection_ulong;
typedef  FArraySection< std::size_t >         FArraySection_size_t;
typedef  FArraySection< std::size_t >         FArraySection_size;
typedef  FArraySection< float >               FArraySection_float;
typedef  FArraySection< double >              FArraySection_double;
typedef  FArraySection< long double >         FArraySection_longdouble;
typedef  FArraySection< char >                FArraySection_char;
typedef  FArraySection< unsigned char >       FArraySection_uchar;
typedef  FArraySection< signed char >         FArraySection_schar;
typedef  FArraySection< std::string >         FArraySection_string;
typedef  FArraySection< Fstring >             FArraySection_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArraySection_fwd_HH
