// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/util.cc
/// @brief
/// @author James Thompson

#include <numeric/kdtree/constants.hh>
#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/construct_kdtree.hh>
#include <numeric/kdtree/KDPoint.hh>
#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/HyperRectangle.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <algorithm>

namespace numeric {
namespace kdtree {

KDNodeOP construct_kd_tree(
	utility::vector1< KDPointOP > & points,

	numeric::Size depth,
	KDTree & tree
) {
	using numeric::Real;
	using numeric::Size;
	using utility::vector1;
	if ( points.size() == 0 ) {
		return nullptr;
	}

	numeric::Size width = points.front()->size();
	numeric::Size const axis( ( depth % width ) + 1 );

	// sort points by this axis, which is the index into the points for
	// comparison. This is a big place where we can get a speedup in tree
	// construction. Optimizations for the future could be:
	// - a linear-time median finding function (smart and hard to
	// implement).
	// - finding a median based on a subset of points and grabbing a pivot
	// point near it (stupid but easy to implement).
	std::sort( points.begin(), points.end(), CompareKDPoints(axis) );

	//std::cout << "axis = " << axis << std::endl;
	//std::cout << "constructing tree from points: " << std::endl;
	//for ( Size ii = 1; ii <= points.size(); ++ii ) {
	// print_point( std::cout, points[ii]->location() );
	// std::cout << std::endl;
	//}
	//std::cout << std::endl << std::endl;

	// location is the median of points along the split axis
	Size const median_idx( static_cast< Size > ( points.size() / 2 ) + 1 );

	KDNodeOP current( new KDNode(tree) );
	//tree.extend_bounds( points[median_idx]->location() );
	//std::cout << "median_idx = " << median_idx << std::endl;
	//std::cout << "current point = ";
	//print_point( std::cout, points[median_idx]->location() );
	//std::cout << std::endl;
	current->point( points[median_idx] );
	current->split_axis( axis );

	// if we have enough points, split the points into two halves:
	// - the set of points less than this location (left_child)
	// - the set of points greater than this location (right_child)
	// Since the points are already sorted by axis.
	if ( points.size() > 1 ) {
		vector1< KDPointOP > left_points(
			points.begin(), points.begin() + median_idx - 1
		);
		current->left_child( construct_kd_tree( left_points, depth + 1, tree ) );
		current->left_child()->parent( current );
	}
	if ( points.size() > 2 ) {
		vector1< KDPointOP > right_points(
			points.begin() + median_idx, points.end()
		);
		current->right_child( construct_kd_tree( right_points, depth + 1, tree ) );
		current->right_child()->parent( current );
	}

	return current;
} // construct_kd_tree

} // kdtree
} // numeric
