// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/kdtree.hh
/// @brief functions for creating a kdtree
/// @author James Thompson

#ifndef INCLUDED_numeric_kdtree_construct_kdtree_hh
#define INCLUDED_numeric_kdtree_construct_kdtree_hh

#include <numeric/types.hh>

#include <numeric/kdtree/KDTree.fwd.hh>
#include <numeric/kdtree/KDPoint.fwd.hh>
#include <numeric/kdtree/KDNode.fwd.hh>

#include <utility/vector1.fwd.hh>

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

} // kdtree
} // numeric

#endif
