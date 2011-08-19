#ifndef INCLUDED_ObjexxFCL_FArray2A_fwd_hh
#define INCLUDED_ObjexxFCL_FArray2A_fwd_hh


// FArray2A Forward Declarations
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
template< typename > class FArray2A;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray2A< bool >                FArray2A_bool;
typedef  FArray2A< byte >                FArray2A_byte;
typedef  FArray2A< sbyte >               FArray2A_sbyte;
typedef  FArray2A< ubyte >               FArray2A_ubyte;
typedef  FArray2A< short int >           FArray2A_short;
typedef  FArray2A< int >                 FArray2A_int;
typedef  FArray2A< long int >            FArray2A_long;
typedef  FArray2A< unsigned short int >  FArray2A_ushort;
typedef  FArray2A< unsigned int >        FArray2A_uint;
typedef  FArray2A< unsigned long int >   FArray2A_ulong;
typedef  FArray2A< std::size_t >         FArray2A_size_t;
typedef  FArray2A< std::size_t >         FArray2A_size;
typedef  FArray2A< float >               FArray2A_float;
typedef  FArray2A< double >              FArray2A_double;
typedef  FArray2A< long double >         FArray2A_longdouble;
typedef  FArray2A< char >                FArray2A_char;
typedef  FArray2A< unsigned char >       FArray2A_uchar;
typedef  FArray2A< signed char >         FArray2A_schar;
typedef  FArray2A< std::string >         FArray2A_string;
typedef  FArray2A< Fstring >             FArray2A_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray2A_fwd_HH
