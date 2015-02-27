// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzTransform.fwd.hh
/// @brief  numeric::xyzTransform forward declarations
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyzTransform_fwd_hh
#define INCLUDED_numeric_xyzTransform_fwd_hh

#include <numeric/types.hh>

// C++ headers
#include <cstddef>


namespace numeric {


// Forward
template< typename > class xyzTransform;


// Types
typedef  xyzTransform< float >  Xformf;
typedef  xyzTransform< double > Xform;

typedef  xyzTransform< numeric::Real > xyzTransform_Real;
typedef  xyzTransform< float  > xyzTransform_float;
typedef  xyzTransform< double > xyzTransform_double;

} // namespace numeric


#endif // INCLUDED_numeric_xyzTransform_FWD_HH
