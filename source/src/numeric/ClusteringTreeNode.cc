// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking/ClusteringTreeNode.cc
/// @brief a node of a hierarchical clustering tree
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <numeric/ClusteringTreeNode.hh>

namespace numeric {

/// @details Auto-generated virtual destructor
ClusteringTreeNode::~ClusteringTreeNode() {}

ClusteringTreeNodeOP ClusteringTreeNode::visit_next_leaf() {

	if ( flag_ == 1 ) {
		flag_ = 2;
		return right();
	}
	if ( flag_ == 0 ) {
		flag_ = 1;
		return left();
	}

	return 0;
}

} // numeric

