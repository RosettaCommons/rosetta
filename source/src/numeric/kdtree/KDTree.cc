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


#include <numeric/types.hh>

#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/HyperRectangle.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace kdtree {

/// @details Auto-generated virtual destructor
KDTree::~KDTree() {}

KDTree::KDTree() :
	size_( 0 ),
	root_( /* NULL */ ),
	bounds_( /* NULL */ )
{}

KDTree::KDTree( utility::vector1< utility::vector1< numeric::Real > > & pts ) {
	utility::vector1< KDPointOP > kd_pts( make_points( pts ) );
	//bounds_ = new HyperRectangle( pts[1], pts[1] );
	bounds_ = HyperRectangleOP( new HyperRectangle( pts ) );
	root_   = construct_kd_tree( kd_pts, 1, *this );
	size( pts.size() );
}

KDTree::KDTree(
	utility::vector1< utility::vector1< numeric::Real > > & pts,
	utility::vector1< utility::pointer::ReferenceCountOP > & data
) {
	assert( pts.size() == data.size() );
	utility::vector1< KDPointOP > kd_pts( make_points( pts, data ) );
	//bounds_ = new HyperRectangle( pts[1], pts[1] );
	bounds_ = HyperRectangleOP( new HyperRectangle( pts ) );
	root_   = construct_kd_tree( kd_pts, 1, *this );
	size( pts.size() );
}

numeric::Size KDTree::size() const {
	return size_;
}

numeric::Size KDTree::ndim() const {
	return bounds_->ndim();
}

KDNodeOP KDTree::root() const {
	return root_;
}

HyperRectangleOP KDTree::bounds() const {
	return bounds_;
}

void KDTree::size( numeric::Size axis ) {
	size_ = axis;
}

void KDTree::root( KDNodeOP new_root ) {
	root_ = new_root;
}

void KDTree::extend_bounds( utility::vector1< numeric::Real > const & pt ) {
	//if ( ! bounds() ) {
	// bounds_ = new HyperRectangle( pt, pt );
	//} else {
	// bounds_->extend( pt );
	//}
	bounds_->extend( pt );
}

} // kdtree
} // numeric
