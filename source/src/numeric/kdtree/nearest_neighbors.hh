// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/nearest_neighbors.hh
/// @brief utility functions for kd-tree. See kdtree.hh for more information.
/// @author James Thompson
//
#ifndef INCLUDED_numeric_kdtree_nearest_neighbors_hh
#define INCLUDED_numeric_kdtree_nearest_neighbors_hh

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <numeric/kdtree/KDTree.fwd.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <numeric/kdtree/KDPointList.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

namespace numeric {
namespace kdtree {

/// @brief Returns a KDPointList of the N nearest neighbors from the KDTree to
/// the given input point.
KDPointList
nearest_neighbors(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	Size const wanted
);

KDPointList
nearest_neighbors(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	Size const wanted,
	numeric::Real const max_dist_allowed
);

/// @brief Searches the KDtree for the nearest neigbor to a given input point,
/// returns nearest neighbor and distance-squared to nearest neigbor by
/// reference.
void
nearest_neighbor(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	// returns:
	KDNodeOP & nearest,
	numeric::Real & dist_sq
);

// Recursive functions that do things to kd-trees, not intended for public
// use because they're a little awkward.

/// @brief returns the nearest neighbor to the given point.
/// @details Parameters are (in order):
/// - current: the base of the tree
/// - pt: the point that is being searched against the tree
/// - bounds: hyper-rectangle in k-space that bounds all points in the tree
/// - max_dist_sq: maximum squared distance that we care about.
///
/// - nearest neighbor (returned by reference)
/// - squared distance to the nearest neighbor
void
nearest_neighbor(
	KDNodeOP & current,
	utility::vector1< numeric::Real > const & pt,
	HyperRectangle & bounds,
	numeric::Real max_dist_sq,

	// returns:
	KDNodeOP & nearest,
	numeric::Real & dist_sq
);

/// @brief Recursive function definition for search for a list of the N nearest
/// neighbors, where N is defined as a member variable of the KDPointList
/// object.
void
nearest_neighbors(
	KDNodeOP & current,
	utility::vector1< numeric::Real > const & pt,
	HyperRectangle & bounds,
	numeric::Real max_dist_sq,

	// returns:
	KDPointList & neighbors
);

} // kdtree
} // numeric

#endif
