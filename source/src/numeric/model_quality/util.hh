// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/model_quality/util.hh
/// @brief utilities for calculated model-quality statistics

/// @author James Thompson

#ifndef INCLUDED_numeric_model_quality_util_hh
#define INCLUDED_numeric_model_quality_util_hh

#include <ObjexxFCL/ObjexxFCL.hh>

#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/RmsData.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace model_quality {

/// @brief Function that takes a vector1 of coordinates, calculates their
/// center of mass, and translates the center of mass to (0,0,0).
template< typename T > void
center_atoms_at_origin(
	utility::vector1< xyzVector< T > > & coords
) {

	// find the center of mass
	xyzVector< T > center_of_mass( 0.0, 0.0, 0.0 );
	typedef typename utility::vector1< xyzVector< T > >::iterator iterator;
	typedef typename utility::vector1< xyzVector< T > >::const_iterator const_iterator;
	for ( const_iterator it = coords.begin(), end = coords.end();
				it != end; ++it
	) {
		center_of_mass += *it;
	}

	center_of_mass /= coords.size();

	// translate all coordinates so that center of mass lies at the origin.
	for ( iterator it = coords.begin(), end = coords.end();
				it != end; ++it
	) {
		*it = *it - center_of_mass;
	}
} // center_atoms_at_origin

} // namespace model_quality
} // namespace numeric

#endif
