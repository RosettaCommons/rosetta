// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/kdtree/KDTree.hh
/// @brief Implementation of a tree in a kd-tree. See numeric/kdtree/kdtree.hh
/// for more information.
/// @author James Thompson
//
#ifndef INCLUDED_numeric_kdtree_KDTree_hh
#define INCLUDED_numeric_kdtree_KDTree_hh

#include <numeric/types.hh>

#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/KDTree.fwd.hh>
#include <numeric/kdtree/HyperRectangle.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

namespace numeric {
namespace kdtree {

class KDTree : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~KDTree();


	/// @brief Empty constructor.
	KDTree();

	/// @brief Constructs a balanced kd-tree from the set of k-dimensional
	/// input points.
	KDTree(
		utility::vector1< utility::vector1< numeric::Real > > & pts
	);

	KDTree(
		utility::vector1< utility::vector1< numeric::Real > > & pts,
		utility::vector1< utility::pointer::ReferenceCountOP > & data
	);

	/// @brief Number of points in the kd-tree.
	numeric::Size size() const;

	/// @brief Number of dimensions in the kd-tree. This is the "k" in kd.
	numeric::Size ndim() const;

	/// @brief Returns the KDNodeOP that is the root of the balanced kd-tree.
	KDNodeOP root() const;

	/// @brief Returns the HyperRectangle that bounds all of the points in the
	/// kd-tree.
	/// @details A HyperRectangle is defined as two vectors upper and lower, with
	/// each dimension of lower having the minimum value seen in each dimension,
	/// and each dimension of higher having the maximum value seen in each
	/// dimension.
	HyperRectangleOP bounds() const;

	/// @brief Sets the number of points in this kd-tree.
	void size( numeric::Size new_size );

	/// @brief Sets the root of the kd-tree.
	void root( KDNodeOP new_root );

	/// @brief Pushes out the bounds of the HyperRectangle bounding this kd-tree
	/// if necessary.
	void extend_bounds( utility::vector1< numeric::Real > const & pt );

private:
	numeric::Size size_;
	KDNodeOP root_;
	HyperRectangleOP bounds_;
}; // class KDTree

} // kdtree
} // numeric

#endif
