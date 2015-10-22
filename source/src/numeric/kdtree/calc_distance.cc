// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/kdtree/calc_distance.cc
/// @brief
/// @author James Thompson

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <cmath>

namespace numeric {
namespace kdtree {

numeric::Real sq_vec_distance(
	utility::vector1< numeric::Real > const & vec1,
	utility::vector1< numeric::Real > const & vec2
) {
	assert( vec1.size() == vec2.size() );
	using numeric::Real;
	using numeric::Size;
	using utility::vector1;

	static Size n_comp( 0 );

	numeric::Real dist( 0.0 );
	for ( vector1< Real >::const_iterator
			it1 = vec1.begin(), it2 = vec2.begin(),
			end1 = vec1.end(), end2 = vec2.end();
			it1 != end1 && it2 != end2; ++it1, ++it2
			) {
		//std::cout << "comparing " << *it1 << " and " << *it2 << std::endl;
		dist += ( *it1 - *it2 ) * ( *it1 - *it2 );
	}
	//std::cout << "dist = " << dist << std::endl;

	n_comp++;
	//std::cout << "n_comp = " << n_comp << std::endl;

	return dist;
}

numeric::Real vec_distance(
	utility::vector1< numeric::Real > const & vec1,
	utility::vector1< numeric::Real > const & vec2
) {
	return std::sqrt( sq_vec_distance( vec1, vec2 ) );
}

} // kdtree
} // numeric
