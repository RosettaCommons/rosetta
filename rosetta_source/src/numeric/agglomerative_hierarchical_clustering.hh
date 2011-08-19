// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/agglomerative_hierarchical_clustering.hh
/// @brief hierarchical clustering utilities
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_numeric_agglomerative_hierarchical_clustering_hh
#define INCLUDED_numeric_agglomerative_hierarchical_clustering_hh

#include <numeric/ClusteringTreeNode.hh>
#include <utility/vector1.hh>

namespace numeric {

template <class T>
void get_cluster_data(utility::vector1<T> &data_in,ClusteringTreeNodeOP cluster,utility::vector1<T> & data_out) {

    utility::vector1<Size> ids;
    cluster->copy_leaf_ids( ids );
    for(Size i=1;i<=ids.size();i++) {
	data_out.push_back( data_in[ids[i]] );
    }
}

utility::vector1<ClusteringTreeNodeOP> single_link_clustering(utility::vector1< utility::vector1<Real> > &,numeric::Size);

utility::vector1<ClusteringTreeNodeOP> average_link_clustering(utility::vector1< utility::vector1<Real> > &,numeric::Size);

utility::vector1<ClusteringTreeNodeOP> complete_link_clustering(utility::vector1< utility::vector1<Real> > &,numeric::Size);

} // numeric

#endif


