// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzTriple.fwd.hh
/// @brief  numeric::xyzTriple forward declarations
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyzTriple_fwd_hh
#define INCLUDED_numeric_xyzTriple_fwd_hh


// C++ headers
#include <cstddef>


namespace numeric {


// Forward
template< typename > class xyzTriple;


// Types
typedef  xyzTriple< bool >                xyzTriple_bool;
typedef  xyzTriple< short int >           xyzTriple_short;
typedef  xyzTriple< int >                 xyzTriple_int;
typedef  xyzTriple< long int >            xyzTriple_long;
typedef  xyzTriple< unsigned short int >  xyzTriple_ushort;
typedef  xyzTriple< unsigned int >        xyzTriple_uint;
typedef  xyzTriple< unsigned long int >   xyzTriple_ulong;
typedef  xyzTriple< std::size_t >         xyzTriple_size_t;
typedef  xyzTriple< std::size_t >         xyzTriple_size;
typedef  xyzTriple< float >               xyzTriple_float;
typedef  xyzTriple< double >              xyzTriple_double;
typedef  xyzTriple< long double >         xyzTriple_longdouble;
typedef  xyzTriple< char >                xyzTriple_char;
typedef  xyzTriple< unsigned char >       xyzTriple_uchar;
typedef  xyzTriple< signed char >         xyzTriple_schar;


} // namespace numeric


#endif // INCLUDED_numeric_xyzTriple_FWD_HH
