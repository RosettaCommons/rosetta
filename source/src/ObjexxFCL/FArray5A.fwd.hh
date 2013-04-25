#ifndef INCLUDED_ObjexxFCL_FArray5A_fwd_hh
#define INCLUDED_ObjexxFCL_FArray5A_fwd_hh


// FArray5A Forward Declarations
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
template< typename > class FArray5A;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray5A< bool >                FArray5A_bool;
typedef  FArray5A< byte >                FArray5A_byte;
typedef  FArray5A< sbyte >               FArray5A_sbyte;
typedef  FArray5A< ubyte >               FArray5A_ubyte;
typedef  FArray5A< short int >           FArray5A_short;
typedef  FArray5A< int >                 FArray5A_int;
typedef  FArray5A< long int >            FArray5A_long;
typedef  FArray5A< unsigned short int >  FArray5A_ushort;
typedef  FArray5A< unsigned int >        FArray5A_uint;
typedef  FArray5A< unsigned long int >   FArray5A_ulong;
typedef  FArray5A< std::size_t >         FArray5A_size_t;
typedef  FArray5A< std::size_t >         FArray5A_size;
typedef  FArray5A< float >               FArray5A_float;
typedef  FArray5A< double >              FArray5A_double;
typedef  FArray5A< long double >         FArray5A_longdouble;
typedef  FArray5A< char >                FArray5A_char;
typedef  FArray5A< unsigned char >       FArray5A_uchar;
typedef  FArray5A< signed char >         FArray5A_schar;
typedef  FArray5A< std::string >         FArray5A_string;
typedef  FArray5A< Fstring >             FArray5A_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray5A_fwd_HH
