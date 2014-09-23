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

private:
    ClusteringTreeNode() {}

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ClusteringTreeNode();

    /// @brief Creates a node with no leaves
    /// @detailed leaves are NULLs, parent is set to this
    static ClusteringTreeNodeOP
    newClusteringTreeNode() {
      ClusteringTreeNodeOP new_node( new ClusteringTreeNode() );
      new_node->this_weak_ptr_ = ClusteringTreeNodeAP( new_node );
      return new_node;
    }

    static ClusteringTreeNodeOP
    newClusteringTreeNode(Size id) {
      ClusteringTreeNodeOP new_node = newClusteringTreeNode();
      new_node->id_ = id;
      new_node->parent_ = new_node->this_weak_ptr_;
      new_node->size_ = 1;
      new_node->distance_ = 0.0;
      new_node->left_.reset();
      new_node->right_.reset();
      return new_node;
    }

    /// @brief Creates a node with given leaves
    /// @detailed parent of the newly created node is set to itself (this pointer); left and right nodes are also chilred of this
    static ClusteringTreeNodeOP
    newClusteringTreeNode(Size id,ClusteringTreeNodeOP left,ClusteringTreeNodeOP right,Real distance = 0.0) {
			ClusteringTreeNodeOP new_node = newClusteringTreeNode();
      new_node->left_ = left;
      new_node->right_ = right;
      new_node->id_ = id;
      new_node->distance_ = distance;
      new_node->size_ = left->size_ + right->size_;
      left->parent_ = new_node->this_weak_ptr_;
      right->parent_ = new_node->this_weak_ptr_;
//      std::cerr<< "Cluster "<<id_<<" created from "<<left_->id_<<" and "<<right_->id_<<" dist: "<<distance<<", size: "<<size_<<"\n";
      return new_node;
    }

    void reset_all_flags() {

      flag_ = 0;

		  ClusteringTreeNodeOP left_op = left(), right_op = right();
      if(left_op) left_op->reset_all_flags();
      if(right_op) right_op->reset_all_flags();
    }

    void set_all_flags(Size new_flag_value) {

		  ClusteringTreeNodeOP left_op = left(), right_op = right();
      if(left_op) left_op->set_all_flags(new_flag_value);
      if(right_op) right_op->set_all_flags(new_flag_value);
      flag_ = new_flag_value;
    }

    ClusteringTreeNodeOP left() { return left_; }

    ClusteringTreeNodeOP right() { return right_; }

    ClusteringTreeNodeOP parent() { return parent_.lock(); }

    Real distance() { return distance_; }

    Size size() { return size_; }

    Size id() { return id_; }

    void left(ClusteringTreeNodeAP new_left) { left_ = new_left.lock(); }

    void right(ClusteringTreeNodeAP new_right) { left_ = new_right.lock(); }

    void parent(ClusteringTreeNodeAP new_parent) { parent_ = new_parent; }

    bool was_visited() { return flag_ == 2; }

    void set_visited() { flag_ = 2; }

    ClusteringTreeNodeOP visit_next_leaf();

    void copy_member_ids(utility::vector1<Size> & dst) {

      dst.push_back(id_);
		  ClusteringTreeNodeOP left_op = left(), right_op = right();
      if(left_op) left_op->copy_member_ids(dst);
      if(right_op) right_op->copy_member_ids(dst);
    }

    void copy_leaf_ids(utility::vector1<Size> & dst) {

      if(left_ || right_) {
        ClusteringTreeNodeOP left_op = left(), right_op = right();
        if(left_op) left_op->copy_leaf_ids(dst);
        if(right_op) right_op->copy_leaf_ids(dst);
      } else {
        dst.push_back(id_);
      }
    }

private:
    ClusteringTreeNodeAP this_weak_ptr_;
    ClusteringTreeNodeAP parent_;
    ClusteringTreeNodeOP left_;
    ClusteringTreeNodeOP right_;
    Size flag_; /// 0 - not visited; 1 - visited 1; 2 - visited 1 & 2
    Size size_;
    Size id_;
    Real distance_;
};

} // numeric

#endif // INCLUDED_numeric_ClusteringTreeNode_HH


