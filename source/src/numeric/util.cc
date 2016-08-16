// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/util.cc
/// @brief small bundle of utilities for dealing with numbers.
/// @author James Thompson

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <algorithm>

namespace numeric {

numeric::Real median( utility::vector1< numeric::Real > const & values ) {
	assert( values.size() ); // An empty list doesn't have a median
	utility::vector1< numeric::Real > vals = values;
	std::sort( vals.begin(), vals.end() );

	numeric::Size const n_vals( vals.size() );
	numeric::Real retval( 0.0 );
	if ( n_vals % 2 == 0 ) { // Even number of items
		retval += 0.5 * vals[ n_vals / 2 ];
		retval += 0.5 * vals[ n_vals / 2 + 1];
	} else { // Odd number of items
		retval = vals[ (n_vals - 1) / 2 + 1 ];
	}
	return retval;
}

numeric::Real mean( utility::vector1< numeric::Real > const & values ) {
	typedef utility::vector1< numeric::Real >::const_iterator iter;

	numeric::Size const n_vals( values.size() );
	numeric::Real total( 0.0 );
	for ( iter it = values.begin(), end = values.end(); it != end; ++it ) {
		total += *it;
	}

	return static_cast< numeric::Real > ( total / n_vals );
}

} // numeric
