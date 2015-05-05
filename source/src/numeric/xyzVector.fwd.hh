// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzVector.fwd.hh
/// @brief  numeric::xyzVector forward declarations
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyzVector_fwd_hh
#define INCLUDED_numeric_xyzVector_fwd_hh


// C++ headers
#include <cstddef>


namespace numeric {


// Forward
template< typename > class xyzVector;


// Types
typedef  xyzVector< bool >                xyzVector_bool;
typedef  xyzVector< short int >           xyzVector_short;
typedef  xyzVector< int >                 xyzVector_int;
typedef  xyzVector< long int >            xyzVector_long;
typedef  xyzVector< unsigned short int >  xyzVector_ushort;
typedef  xyzVector< unsigned int >        xyzVector_uint;
typedef  xyzVector< unsigned long int >   xyzVector_ulong;
typedef  xyzVector< std::size_t >         xyzVector_size_t;
typedef  xyzVector< std::size_t >         xyzVector_size;
typedef  xyzVector< float >               xyzVector_float;
typedef  xyzVector< double >              xyzVector_double;
typedef  xyzVector< long double >         xyzVector_longdouble;
typedef  xyzVector< char >                xyzVector_char;
typedef  xyzVector< unsigned char >       xyzVector_uchar;
typedef  xyzVector< signed char >         xyzVector_schar;


} // namespace numeric

#endif // INCLUDED_numeric_xyzVector_FWD_HH
