#ifndef INCLUDED_ObjexxFCL_ChunkVector_fwd_hh
#define INCLUDED_ObjexxFCL_ChunkVector_fwd_hh


// ChunkVector Forward Declarations
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
template< typename > class ChunkVector;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  ChunkVector< bool >                ChunkVector_bool;
typedef  ChunkVector< byte >                ChunkVector_byte;
typedef  ChunkVector< sbyte >               ChunkVector_sbyte;
typedef  ChunkVector< ubyte >               ChunkVector_ubyte;
typedef  ChunkVector< short int >           ChunkVector_short;
typedef  ChunkVector< int >                 ChunkVector_int;
typedef  ChunkVector< long int >            ChunkVector_long;
typedef  ChunkVector< unsigned short int >  ChunkVector_ushort;
typedef  ChunkVector< unsigned int >        ChunkVector_uint;
typedef  ChunkVector< unsigned long int >   ChunkVector_ulong;
typedef  ChunkVector< std::size_t >         ChunkVector_size_t;
typedef  ChunkVector< std::size_t >         ChunkVector_size;
typedef  ChunkVector< float >               ChunkVector_float;
typedef  ChunkVector< double >              ChunkVector_double;
typedef  ChunkVector< long double >         ChunkVector_longdouble;
typedef  ChunkVector< char >                ChunkVector_char;
typedef  ChunkVector< unsigned char >       ChunkVector_uchar;
typedef  ChunkVector< signed char >         ChunkVector_schar;
typedef  ChunkVector< std::string >         ChunkVector_string;
typedef  ChunkVector< Fstring >             ChunkVector_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ChunkVector_fwd_HH
