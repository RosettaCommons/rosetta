// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/kdtree/KDNode.hh
/// @brief Implementation of a node in a kd-tree. See numeric/kdtree/kdtree.hh
/// for more information.
/// @author James Thompson
//
#ifndef INCLUDED_numeric_kdtree_KDNode_hh
#define INCLUDED_numeric_kdtree_KDNode_hh

#include <numeric/types.hh>

#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/KDTree.fwd.hh>
#include <numeric/kdtree/KDPoint.fwd.hh>
#include <numeric/kdtree/KDPoint.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace kdtree {

class KDNode : public utility::pointer::ReferenceCount {

public:
	/// @brief Constructor for a KDNode. Takes a const-refence
	/// a KDTree that should contain this KDNode.
	KDNode( KDTree const & tree );
	virtual ~KDNode();

	/// @brief Returns the parent of this KDNode in the tree,
	/// or NULL if there is no parent.
	KDNodeOP parent() const;

	/// @brief Returns the left child of this KDNode in the tree,
	/// or NULL if there is no left child.
	KDNodeOP left_child() const;

	/// @brief Returns the right child of this KDNode in the tree,
	/// or NULL if there is no right child.
	KDNodeOP right_child() const;

	/// @brief Returns a const reference to the Tree that contains
	/// this KDNode.
	KDTree const & tree() const;

	/// @brief Returns the location of this node in k-space.
	utility::vector1< numeric::Real > location() const;

	utility::pointer::ReferenceCountOP data() const;

	/// @brief Returns the dimension along which this node splits points.
	numeric::Size split_axis() const;

	/// @brief Returns true if this node has no children, false otherwise.
	bool is_leaf() const;

	/// @brief Returns true if this node has no parent, false otherwise.
	bool is_root() const;

	KDPointOP point() const;

	numeric::Real distance() const;

	/// @brief Sets the parent for this node.
	void parent( KDNodeOP new_parent );

	/// @brief Sets the left child for this node.
	void left_child( KDNodeOP new_left_child );

	/// @brief Sets the right child for this node.
	void right_child( KDNodeOP new_right_child );

	/// @brief Sets the location of this node in k-space.
	void location( utility::vector1< numeric::Real > new_location );

	void point( KDPointOP new_point );

	void distance( numeric::Real new_dist );

	/// @brief Sets the dimension along which this node splits points.
	void split_axis( numeric::Size axis );

	// Undefined, commenting out
	/// @brief output operator for KDNode
	/*
	friend std::ostream & operator<< (
	std::ostream & out,
	const KDNode & kdnode
	); */

private:
	numeric::Size split_axis_;
	KDNodeOP parent_, left_child_, right_child_;
	KDPointOP pt_;
	KDTree const & tree_;
}; // class KDNode


} // kdtree
} // numeric

#endif
