// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/kdtree.hh
/// @brief utility functions for kd-tree. See kdtree.hh for more information.
/// @author James Thompson
#ifndef INCLUDED_numeric_kdtree_calc_distance_hh
#define INCLUDED_numeric_kdtree_calc_distance_hh

#include <numeric/types.hh>
#include <utility/vector1.fwd.hh>

namespace numeric {
namespace kdtree {

/// distance metrics for real-valued points

/// @brief Returns the square of the Euclidean distance between the two points
/// vec1 and vec2.
numeric::Real sq_vec_distance(
	utility::vector1< numeric::Real > const & vec1,
	utility::vector1< numeric::Real > const & vec2
);

/// @brief Returns the Euclidean distance between the two points vec1 and vec2.
numeric::Real vec_distance(
	utility::vector1< numeric::Real > const & vec1,
	utility::vector1< numeric::Real > const & vec2
);

} // kdtree
} // numeric

#endif
