// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  numeric/ClusteringTreeNode.hh
/// @brief a node of a hierarchical clustering tree
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_numeric_ClusteringTreeNode_hh
#define INCLUDED_numeric_ClusteringTreeNode_hh

#include <numeric/ClusteringTreeNode.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/types.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>

namespace numeric {


class ClusteringTreeNode : public utility::pointer::ReferenceCount {

public:

    /// @brief Creates a node with no leaves
    /// @detailed leaves are NULLs, parent is set to this
    ClusteringTreeNode(Size id) {

	id_ = id;
	parent_ = this;
	size_ = 1;
	distance_ = 0.0;
	left_ = right_ = 0;
    }

    /// @brief Creates a node with given leaves
    /// @detailed parent of the newly created node is set to itself (this pointer); left and right nodes are also chilred of this
    ClusteringTreeNode(Size id,ClusteringTreeNodeOP left,ClusteringTreeNodeOP right,Real distance = 0.0) {

	parent_ = this;
	left_ = left;
	left_->parent_ = this;
	right_ = right;
	right_->parent_ = this;
	id_ = id;
	distance_ = distance;
	size_ = left->size_ + right_->size_;
//	std::cerr<< "Cluster "<<id_<<" created from "<<left_->id_<<" and "<<right_->id_<<" dist: "<<distance<<", size: "<<size_<<"\n";
    }

    void reset_all_flags() {

	flag_ = 0;
	if(left_) left_->reset_all_flags();
	if(right_) right_->reset_all_flags();
    }

    void set_all_flags(Size new_flag_value) {

	if(left_) left_->set_all_flags(new_flag_value);
	if(right_) right_->set_all_flags(new_flag_value);
	flag_ = new_flag_value;
    }

    ClusteringTreeNodeOP left() { return left_; }

    ClusteringTreeNodeOP right() { return right_; }

    ClusteringTreeNodeOP parent() { return parent_; }

    Real distance() { return distance_; }

    Size size() { return size_; }

    Size id() { return id_; }

    void left(ClusteringTreeNodeOP new_left) { left_ = new_left; }

    void right(ClusteringTreeNodeOP new_right) { left_ = new_right; }

    void parent(ClusteringTreeNodeOP new_parent) { parent_ = new_parent; }

    bool was_visited() { return flag_ == 2; }

    void set_visited() { flag_ = 2; }

    ClusteringTreeNodeOP visit_next_leaf();

    void copy_member_ids(utility::vector1<Size> & dst) {

	dst.push_back(id_);
	if(left_) left_->copy_member_ids(dst);
	if(right_) right_->copy_member_ids(dst);
    }

    void copy_leaf_ids(utility::vector1<Size> & dst) {

	if(left_ || right_) {
	    if(left_) left_->copy_leaf_ids(dst);
	    if(right_) right_->copy_leaf_ids(dst);
	} else {
	    dst.push_back(id_);
	}
    }

private:
    ClusteringTreeNodeOP parent_;
    ClusteringTreeNodeOP left_;
    ClusteringTreeNodeOP right_;
    Size flag_; /// 0 - not visited; 1 - visited 1; 2 - visited 1 & 2
    Size size_;
    Size id_;
    Real distance_;
};

} // numeric

#endif // INCLUDED_numeric_ClusteringTreeNode_HH


