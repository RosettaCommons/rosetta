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
#include <numeric/types.hh>

namespace numeric {

class ClusterOptions {
	public:
	
	ClusterOptions() :
	node1_( 0 ),
	node2_( 0 ),
	min_i_( 0 ),
	min_j_( 0 ),
	min_  ( 0 ) {}
	
	ClusterOptions( Size node1, Size node2, Size min_i, Size min_j, Real min ) :
	node1_( node1 ),
	node2_( node2 ),
	min_i_( min_i ),
	min_j_( min_j ),
	min_  ( min   ) {}
	
	Size node1_;
	Size node2_;
	Size min_i_;
	Size min_j_;
	Real min_;
};

template <class T>
void get_cluster_data(utility::vector1<T> &data_in,ClusteringTreeNodeOP cluster,utility::vector1<T> & data_out) {

	utility::vector1<Size> ids;
	cluster->copy_leaf_ids( ids );
	for(Size i=1;i<=ids.size();i++) {
		data_out.push_back( data_in[ids[i]] );
	}
}


class AgglomerativeHierarchicalClusterer {
public:
	utility::vector1<ClusteringTreeNodeOP>cluster( utility::vector1< utility::vector1<Real> > & dm, Size n );
	virtual void comparator(utility::vector1< utility::vector1<Real> > &distance_matrix,
							utility::vector1<Size> const & members1,
							utility::vector1<Size> const & members2,
							ClusterOptions & co ) = 0; /*Size node1, Size node2,
							Size& min_i, Size& min_j,
							Real& min ) = 0;*/
};

class SingleLinkClusterer : public AgglomerativeHierarchicalClusterer {
public:
	void comparator( utility::vector1< utility::vector1<Real> > &distance_matrix,
					utility::vector1<Size> const & members1,
					utility::vector1<Size> const & members2,
					ClusterOptions & co);
};

class AverageLinkClusterer : public AgglomerativeHierarchicalClusterer {
public:
	void comparator( utility::vector1< utility::vector1<Real> > &distance_matrix,
					utility::vector1<Size> const & members1,
					utility::vector1<Size> const & members2,
					ClusterOptions & co);
};

class CompleteLinkClusterer : public AgglomerativeHierarchicalClusterer {
public:
	void comparator( utility::vector1< utility::vector1<Real> > &distance_matrix,
					utility::vector1<Size> const & members1,
					utility::vector1<Size> const & members2,
					ClusterOptions & co);
};

} // numeric

#endif
