// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/deriv/angle_deriv.hh
/// @brief  inline function for computing f1/f2 derivatives for a function of a distance
/// @author Phil Bradley did all the hard work deriving the math represented here.
/// @author Andrew Leaver-Fay copy-and-pasted Phil's code into this file from
/// the AtomPairConstraint.cc file for general use.

#ifndef INCLUDED_numeric_deriv_distance_deriv_hh
#define INCLUDED_numeric_deriv_distance_deriv_hh

#include <numeric/xyzVector.hh>

namespace numeric {
namespace deriv   {

/// @brief Compute the f1/f2 derivative vectors for point p1 for a function F of the
/// distance between p1 and p2.  This function returns the distance which should
/// be used to evaluate dF_ddist.  dF_ddist should then be multiplied into both
/// f1 and f2.  The values of the output variables f1 and f2 are overwritten.
template < class P >
inline
void
distance_f1_f2_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	P & distance,
	xyzVector< P > & f1,
	xyzVector< P > & f2
)
{
	typedef P Real;

	f2 = p1 - p2;
	distance = f2.length();
	if ( distance != Real(0.0) ) {
		Real const invd = Real(1.0) / distance;
		f1 = p1.cross( p2 );
		f1 *= invd;
		f2 *= invd;
	} else {
		f1 = Real(0.0);
	}
}


}
}

#endif
