// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/numeric:agglomerative_hierarchical_clustering
/// @brief hierarchical clustering utilities
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <numeric/ClusteringTreeNode.hh>

//#include <basic/Tracer.hh>

#include <list>
#include <iostream>

namespace numeric {

utility::vector1<ClusteringTreeNodeOP>
    single_link_clustering(utility::vector1< utility::vector1<Real> > &distance_matrix,Size n_clusters) {

    std::list<Size> active_indexes;
    utility::vector1<ClusteringTreeNodeOP> nodes;
    int cluster_id = 1;
    ClusteringTreeNodeOP root;

    for(Size i=1;i<=distance_matrix.size();i++) {
	active_indexes.push_back(i);
	nodes.push_back( ClusteringTreeNode::newClusteringTreeNode(i) );
	cluster_id++;
    }

    while(active_indexes.size() >n_clusters) {
        active_indexes.sort();
        Size min_i = 0;
	Size min_j = 0;
	Real min = 99999999999.0;
	std::list<Size>::iterator it1,it2;

	for(std::list<Size>::iterator it1=active_indexes.begin() ; it1 != active_indexes.end(); it1++) {
	    utility::vector1<Size> members1;
	    Size node1 = *it1;
            nodes[node1]->copy_leaf_ids( members1 );
	    for(std::list<Size>::iterator it2=active_indexes.begin() ; it2 != active_indexes.end(); it2++) {
	        Size node2 = *it2;
	        if(node1<=node2)
	    	    break;
    		utility::vector1<Size> members2;
        	nodes[node2]->copy_leaf_ids( members2 );
	        for(utility::vector1<Size>::iterator m1=members1.begin(); m1 !=members1.end();m1++) {
    		    for(utility::vector1<Size>::iterator m2=members2.begin(); m2 !=members2.end();m2++) {
			if(distance_matrix[*m1][*m2] < min) {
			    min = distance_matrix[*m1][*m2];
			    min_i = node1;
			    min_j = node2;
			}
		    }
		}
	    }
	}
	root = ClusteringTreeNode::newClusteringTreeNode(cluster_id,nodes[min_i],nodes[min_j],min);
	cluster_id++;
	nodes.push_back(root);
	active_indexes.remove(min_i);
	active_indexes.remove(min_j);
	active_indexes.push_back(nodes.size());
//	std::cerr<<active_indexes.size() <<" "<<nodes.size()<<"\n";
    }

    assert(active_indexes.size() == n_clusters);
    utility::vector1<ClusteringTreeNodeOP> out;
    for(std::list<Size>::iterator it=active_indexes.begin() ; it != active_indexes.end(); it++)
	out.push_back(nodes[*it]);

    return out;
}


utility::vector1<ClusteringTreeNodeOP>
    average_link_clustering(utility::vector1< utility::vector1<Real> > &distance_matrix,Size n_clusters) {

    std::list<Size> active_indexes;
    utility::vector1<ClusteringTreeNodeOP> nodes;
    int cluster_id = 1;
    ClusteringTreeNodeOP root;

    for(Size i=1;i<=distance_matrix.size();i++) {
	active_indexes.push_back(i);
	nodes.push_back( ClusteringTreeNode::newClusteringTreeNode(i) );
	cluster_id++;
    }

    while(active_indexes.size() >n_clusters) {
        active_indexes.sort();
        Size min_i = 0;
	Size min_j = 0;
	Real min = 99999999999.0;
	std::list<Size>::iterator it1,it2;

	for(std::list<Size>::iterator it1=active_indexes.begin() ; it1 != active_indexes.end(); it1++) {
	    utility::vector1<Size> members1;
	    Size node1 = *it1;
            nodes[node1]->copy_leaf_ids( members1 );
	    for(std::list<Size>::iterator it2=active_indexes.begin() ; it2 != active_indexes.end(); it2++) {
	        Size node2 = *it2;
	        if(node1<=node2)
	    	    break;
    		utility::vector1<Size> members2;
        	nodes[node2]->copy_leaf_ids( members2 );
        	Real dist = 0.0;
        	Real n = 0.0;
	        for(utility::vector1<Size>::iterator m1=members1.begin(); m1 !=members1.end();m1++) {
    		    for(utility::vector1<Size>::iterator m2=members2.begin(); m2 !=members2.end();m2++) {
    			dist += distance_matrix[*m1][*m2];
    			n++;
		    }
		}
		dist = dist / n;
		if(dist < min) {
		    min = dist;
		    min_i = node1;
		    min_j = node2;
		}
	    }
	}
	root = ClusteringTreeNode::newClusteringTreeNode(cluster_id,nodes[min_i],nodes[min_j],min);
	cluster_id++;
	nodes.push_back(root);
	active_indexes.remove(min_i);
	active_indexes.remove(min_j);
	active_indexes.push_back(nodes.size());
//	std::cerr<<active_indexes.size() <<" "<<nodes.size()<<"\n";
    }

    assert(active_indexes.size() == n_clusters);
    utility::vector1<ClusteringTreeNodeOP> out;
    for(std::list<Size>::iterator it=active_indexes.begin() ; it != active_indexes.end(); it++)
	out.push_back(nodes[*it]);

    return out;
}

utility::vector1<ClusteringTreeNodeOP>
    complete_link_clustering(utility::vector1< utility::vector1<Real> > &distance_matrix,Size n_clusters) {

    std::list<Size> active_indexes;
    utility::vector1<ClusteringTreeNodeOP> nodes;
    int cluster_id = 1;
    ClusteringTreeNodeOP root;

    for(Size i=1;i<=distance_matrix.size();i++) {
	active_indexes.push_back(i);
	nodes.push_back( ClusteringTreeNode::newClusteringTreeNode(i) );
	cluster_id++;
    }

    while(active_indexes.size() >n_clusters) {
        active_indexes.sort();
        Size min_i = 0;
	Size min_j = 0;
	Real min = 99999999999.0;
	std::list<Size>::iterator it1,it2;

	for(std::list<Size>::iterator it1=active_indexes.begin() ; it1 != active_indexes.end(); it1++) {
	    utility::vector1<Size> members1;
	    Size node1 = *it1;
            nodes[node1]->copy_leaf_ids( members1 );
	    for(std::list<Size>::iterator it2=active_indexes.begin() ; it2 != active_indexes.end(); it2++) {
	        Size node2 = *it2;
	        if(node1<=node2)
	    	    break;
    		utility::vector1<Size> members2;
        	nodes[node2]->copy_leaf_ids( members2 );
	        Real max = -9999999999.0;
	        for(utility::vector1<Size>::iterator m1=members1.begin(); m1 !=members1.end();m1++) {
    		    for(utility::vector1<Size>::iterator m2=members2.begin(); m2 !=members2.end();m2++) {
			if(distance_matrix[*m1][*m2] > max) {
			    max = distance_matrix[*m1][*m2];
			}
		    }
		}
		if(max <  min) {
		    min = max;
		    min_i = node1;
		    min_j = node2;
		}
	    }
	}
	root = ClusteringTreeNode::newClusteringTreeNode(cluster_id,nodes[min_i],nodes[min_j],min);
	cluster_id++;
	nodes.push_back(root);
	active_indexes.remove(min_i);
	active_indexes.remove(min_j);
	active_indexes.push_back(nodes.size());
    }

    assert(active_indexes.size() == n_clusters);
    utility::vector1<ClusteringTreeNodeOP> out;
    for(std::list<Size>::iterator it=active_indexes.begin() ; it != active_indexes.end(); it++)
	out.push_back(nodes[*it]);

    return out;
}


} // numeric



