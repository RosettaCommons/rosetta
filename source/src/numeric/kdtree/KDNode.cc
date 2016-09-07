// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/KDNode.cc
/// @brief Implementation of a node in a kd-tree. See numeric/kdtree/kdtree.hh
/// for more information.
/// @author James Thompson


#include <numeric/types.hh>

#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/KDPoint.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace kdtree {

KDNode::KDNode( KDTree const & tree ) :
	split_axis_( 0 ),
	parent_( /* NULL */ ),
	left_child_( /* NULL */ ),
	right_child_( /* NULL */ ),
	pt_( /* NULL */ ),
	tree_( tree )
{}

KDNode::~KDNode() = default;

KDNodeOP KDNode::parent() const {
	return parent_;
}

KDNodeOP KDNode::left_child() const {
	return left_child_;
}

KDNodeOP KDNode::right_child() const {
	return right_child_;
}

KDTree const & KDNode::tree() const {
	return tree_;
}

utility::vector1< numeric::Real > KDNode::location() const {
	return pt_->location();
}

utility::pointer::ReferenceCountOP KDNode::data() const {
	return pt_->data();
}

numeric::Size KDNode::split_axis() const {
	return split_axis_;
}

bool KDNode::is_leaf() const {
	return ( !right_child() && !left_child() );
}

bool KDNode::is_root() const {
	return( !parent() );
}

KDPointOP KDNode::point() const {
	return pt_;
}

numeric::Real KDNode::distance() const {
	return pt_->distance();
}

void KDNode::distance( numeric::Real new_dist ) {
	pt_->distance( new_dist );
}

void KDNode::parent( KDNodeOP new_parent ) {
	parent_ = new_parent;
}

void KDNode::left_child( KDNodeOP new_left_child ) {
	left_child_ = new_left_child;
}

void KDNode::right_child( KDNodeOP new_right_child ) {
	right_child_ = new_right_child;
}

void KDNode::point( KDPointOP new_point ) {
	pt_ = new_point;
}

void KDNode::location( utility::vector1< numeric::Real > loc ) {
	pt_->location( loc );
}

void KDNode::split_axis( numeric::Size axis ) {
	split_axis_ = axis;
}

std::ostream & operator<< (
	std::ostream & out,
	const utility::vector1< numeric::Real > & points
) {
	using numeric::Real;
	using utility::vector1;
	for (double point : points) {
		out << ' ' << point;
	}
	return out;
}

} // kdtree
} // numeric
