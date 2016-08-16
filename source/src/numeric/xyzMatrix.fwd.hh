// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzMatrix.fwd.hh
/// @brief  numeric::xyzMatrix forward declarations
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyzMatrix_fwd_hh
#define INCLUDED_numeric_xyzMatrix_fwd_hh

#include <platform/types.hh>

// C++ headers
#include <cstddef>


namespace numeric {


// Forward
template< typename > class xyzMatrix;


// Types
typedef  xyzMatrix< bool >                xyzMatrix_bool;
typedef  xyzMatrix< short int >           xyzMatrix_short;
typedef  xyzMatrix< int >                 xyzMatrix_int;
typedef  xyzMatrix< long int >            xyzMatrix_long;
typedef  xyzMatrix< unsigned short int >  xyzMatrix_ushort;
typedef  xyzMatrix< unsigned int >        xyzMatrix_uint;
typedef  xyzMatrix< unsigned long int >   xyzMatrix_ulong;
typedef  xyzMatrix< platform::Size >      xyzMatrix_Size;
typedef  xyzMatrix< platform::Size >      xyzMatrix_size_t;
typedef  xyzMatrix< platform::Size >      xyzMatrix_size;
typedef  xyzMatrix< float >               xyzMatrix_float;
typedef  xyzMatrix< double >              xyzMatrix_double;
typedef  xyzMatrix< long double >         xyzMatrix_longdouble;
typedef  xyzMatrix< char >                xyzMatrix_char;
typedef  xyzMatrix< unsigned char >       xyzMatrix_uchar;
typedef  xyzMatrix< signed char >         xyzMatrix_schar;


} // namespace numeric


#endif // INCLUDED_numeric_xyzMatrix_FWD_HH

