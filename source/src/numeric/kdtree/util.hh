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
//
#ifndef INCLUDED_numeric_kdtree_util_hh
#define INCLUDED_numeric_kdtree_util_hh

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <numeric/kdtree/KDTree.fwd.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/KDPoint.fwd.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <numeric/kdtree/WrappedReal.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <iostream>

namespace numeric {
namespace kdtree {

/// @brief Function for constructing a KDTree. Returns a KDNodeOP that
/// represents the root of the tree. Points need to be sorted as the
/// tree is being constructed, so the reference to the points is non-const.
KDNodeOP construct_kd_tree(
	utility::vector1< KDPointOP > & points,
	numeric::Size depth,
	KDTree & tree
);

/// distance metrics for real-valued points

/// @brief Transforms the list of points given into percentiles using
/// a linear mapping from the input space to percentile-space for each
/// variable.
/// @details For each variable X in row R, replaces X with the quantity
/// ( X - min(R) ) / ( max(R) - min(R) ). Runs in O(N) time.
void transform_percentile(
	utility::vector1< utility::vector1< numeric::Real > > & points
);

void transform_percentile(
	utility::vector1< utility::vector1< numeric::Real > > & points,
	HyperRectangleOP bounds
);

void transform_percentile_single_pt(
	utility::vector1< numeric::Real > & point,
	HyperRectangleOP bounds
);

/// @brief Makes a vector of KDPoints.
utility::vector1< KDPointOP > make_points(
	utility::vector1< utility::vector1< numeric::Real > > const & points
);

void print_points(
	std::ostream & out,
	utility::vector1< utility::vector1< numeric::Real > > const & points
);

void print_point(
	std::ostream & out,
	utility::vector1< numeric::Real > const & point
);

HyperRectangleOP get_percentile_bounds(
	utility::vector1< utility::vector1< numeric::Real > > & points
);

/// @brief Makes a vector1 of KDPoints, associating the nth entry in data
/// with the nth entry in points.
utility::vector1< KDPointOP > make_points(
	utility::vector1< utility::vector1< numeric::Real > > const & points,
	utility::vector1< utility::pointer::ReferenceCountOP > const & data
);

/// @brief output operator for vector1< Real >
std::ostream & operator<< (
	std::ostream & out,
	const utility::vector1< numeric::Real > & points
);

/// @brief returns true if the given hyper-rectangle intersects with the given
/// hypersphere.
bool hr_intersects_hs(
	HyperRectangle hr,
	utility::vector1< numeric::Real > const & pt,
	numeric::Real const r
);

void print_tree(
	std::ostream & out,
	KDNodeOP const & current,
	Size current_depth,
	Size const width = 3
);

bool is_legal_less_than(
	KDNodeOP const & current,
	Size const split_axis,
	Real const split_value
);

bool is_legal_greater_than(
	KDNodeOP const & current,
	Size const split_axis,
	Real const split_value
);

bool is_legal_kdtree_below_node(
	KDNodeOP const & current
);

} // kdtree
} // numeric

#endif
