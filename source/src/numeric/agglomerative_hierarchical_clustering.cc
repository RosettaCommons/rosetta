// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/numeric:agglomerative_hierarchical_clustering
/// @brief hierarchical clustering utilities
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <numeric/ClusteringTreeNode.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>


//#include <basic/Tracer.hh>

#include <list>

namespace numeric {

utility::vector1<ClusteringTreeNodeOP>
AgglomerativeHierarchicalClusterer::cluster(
	utility::vector1< utility::vector1<Real> > &distance_matrix, Size n_clusters ) {

	std::list<Size> active_indexes;
	utility::vector1<ClusteringTreeNodeOP> nodes;
	int cluster_id = 1;
	ClusteringTreeNodeOP root;

	for ( Size i = 1; i <= distance_matrix.size(); i++ ) {
		active_indexes.push_back(i);
		nodes.push_back( ClusteringTreeNode::newClusteringTreeNode(i) );
		cluster_id++;
	}

	while ( active_indexes.size() >n_clusters ) {
		active_indexes.sort();
		ClusterOptions co( 0, 0, 0, 0, 99999999999.0);

		for ( auto it1 = active_indexes.begin() ; it1 != active_indexes.end(); ++it1 ) {
			utility::vector1<Size> members1;
			Size node1 = *it1;
			nodes[node1]->copy_leaf_ids( members1 );
			for ( unsigned long node2 : active_indexes ) {
				if ( node1 <= node2 ) break;
				utility::vector1< Size > members2;
				nodes[ node2 ]->copy_leaf_ids( members2 );
				co.node1_ = node1; co.node2_ = node2;
				comparator( distance_matrix, members1, members2, co );
			}
		}
		root = ClusteringTreeNode::newClusteringTreeNode(cluster_id,nodes[co.min_i_],nodes[co.min_j_],co.min_);
		cluster_id++;
		nodes.push_back(root);
		active_indexes.remove(co.min_i_);
		active_indexes.remove(co.min_j_);
		active_indexes.push_back(nodes.size());
		// std::cerr<<active_indexes.size() <<" "<<nodes.size()<<"\n";
	}

	assert(active_indexes.size() == n_clusters);
	utility::vector1<ClusteringTreeNodeOP> out;
	for ( unsigned long & active_indexe : active_indexes ) {
		out.push_back(nodes[active_indexe]);
	}

	return out;
}

void SingleLinkClusterer::comparator(
	utility::vector1< utility::vector1<Real> > &distance_matrix,
	utility::vector1<Size> const & members1,
	utility::vector1<Size> const & members2,
	ClusterOptions & co
) {
	for ( unsigned long m1 : members1 ) {
		for ( unsigned long m2 : members2 ) {
			if ( distance_matrix[ m1 ][ m2 ] < co.min_ ) {
				co.min_ = distance_matrix[ m1 ][ m2 ];
				co.min_i_ = co.node1_;
				co.min_j_ = co.node2_;
			}
		}
	}
}

void AverageLinkClusterer::comparator(
	utility::vector1< utility::vector1<Real> > &distance_matrix,
	utility::vector1<Size> const & members1,
	utility::vector1<Size> const & members2,
	ClusterOptions & co
) {
	Real dist = 0.0;
	Real n = 0.0;
	for ( unsigned long m1 : members1 ) {
		for ( unsigned long m2 : members2 ) {
			dist += distance_matrix[ m1 ][ m2 ];
			n++;
		}
	}
	dist = dist / n;
	if ( dist < co.min_ ) {
		co.min_ = dist;
		co.min_i_ = co.node1_;
		co.min_j_ = co.node2_;
	}
}

void CompleteLinkClusterer::comparator(
	utility::vector1< utility::vector1<Real> > &distance_matrix,
	utility::vector1<Size> const & members1,
	utility::vector1<Size> const & members2,
	ClusterOptions & co
) {
	Real max = -9999999999.0;
	for ( unsigned long m1 : members1 ) {
		for ( unsigned long m2 : members2 ) {
			if ( distance_matrix[ m1 ][ m2 ] > max ) {
				max = distance_matrix[ m1 ][ m2 ];
			}
		}
	}
	if ( max < co.min_ ) {
		co.min_ = max;
		co.min_i_ = co.node1_;
		co.min_j_ = co.node2_;
	}
}

} // numeric
