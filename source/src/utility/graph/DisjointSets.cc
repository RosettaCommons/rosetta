// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/Graph.hh
/// @brief  generic graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit Headers
#include <utility/graph/DisjointSets.hh>

namespace utility {
namespace graph {

///
/// @brief
/// default constructor, for when number of nodes is not known
///
DisjointSets::DisjointSets() {}

///
/// @brief
/// constructor for when number of nodes is known up front. fastest.
///
DisjointSets::DisjointSets( platform::Size n_nodes ) :
	nodes_( n_nodes )
{
	for ( platform::Size ii = 1; ii <= n_nodes; ++ii ) {
		nodes_[ ii ].parent = ii;
		nodes_[ ii ].rank = 0;
	}
}

///
/// @brief
/// returns the total number of nodes
///
platform::Size DisjointSets::n_nodes() const {
	return nodes_.size();
}

///
/// @brief
/// creates a new set
///
/// @details This implementation uses indices to objects in arrays instead of pointers and so it relies on vector
/// push-back methods (O(N)) instead of list push-back methods (O(1)).  If enough people clamour, I'll go back and make
/// this faster...(apl)
///
void DisjointSets::ds_make_set() {
	nodes_.resize( n_nodes() + 1 );
	nodes_[ n_nodes() ].parent = n_nodes();
	nodes_[ n_nodes() ].rank = 0;
}

///
/// @brief
/// given a node_id, return the representative for that node
///
platform::Size
DisjointSets::ds_find( platform::Size node_id ) const {

	if ( nodes_[ node_id ].parent != node_id ) {
		nodes_[ node_id ].parent = ds_find( nodes_[ node_id ].parent );
	}
	return nodes_[ node_id ].parent;
}

///
/// @brief
/// combine two sets; make it so that two nodes end up in the same set
///
void
DisjointSets::ds_union( platform::Size node1, platform::Size node2 ) {

	platform::Size parent_node1 = ds_find( node1 );
	platform::Size parent_node2 = ds_find( node2 );

	if ( nodes_[ parent_node1 ].rank < nodes_[ parent_node2 ].rank ) {
		nodes_[ parent_node1 ].parent = parent_node2;
		++nodes_[ parent_node2 ].rank;

	} else if ( nodes_[ parent_node1 ].rank > nodes_[ parent_node2 ].rank ) {
		nodes_[ parent_node2 ].parent = parent_node1;
		++nodes_[ parent_node1 ].rank;

	} else if ( parent_node1 != parent_node2 ) {
		nodes_[ parent_node1 ].parent = parent_node2;
		++nodes_[ parent_node2 ].rank;
	}
}

///
/// @brief
/// count the number of disjoint sets. O(N)
///
platform::Size DisjointSets::n_disjoint_sets() const {

	platform::Size n_disjoint( 0 );
	for ( platform::Size ii = 1; ii <= n_nodes(); ++ii ) {
		if ( nodes_[ ii ].parent == ii ) {
			++n_disjoint;
		}
	}
	return n_disjoint;
}

///
/// @brief
/// returns a vector1 containing the size of each disjoint set. O(N)
///
utility::vector1< platform::Size >
DisjointSets::disjoint_set_sizes() const {

	utility::vector1< platform::Size > index_2_ds( n_nodes(), 0 );
	platform::Size n_disjoint( 0 );
	for ( platform::Size ii = 1; ii <= n_nodes(); ++ii ) {
		if ( nodes_[ ii ].parent == ii ) {
			index_2_ds[ ii ] = ++n_disjoint;
		}
	}

	utility::vector1< platform::Size > ds_set_sizes( n_disjoint, 0 );
	for ( platform::Size ii = 1; ii <= n_nodes(); ++ii ) {
		++ds_set_sizes[ index_2_ds[ ds_find( ii ) ] ];
	}

	return ds_set_sizes;
}

///
/// @brief
/// returns a vector1 of the nodes in the set containing the specified node.
///
utility::vector1< platform::Size >
DisjointSets::nodes_in_set( platform::Size node_id ) const {

	utility::vector1< platform::Size > nis;

	platform::Size const root = ds_find( node_id );

	for ( platform::Size i = 1, ie = n_nodes(); i <= ie; ++i ) {
		if ( root == ds_find( i ) ) {
			nis.push_back( i );
		}
	}

	return nis;
}

///
/// @brief
/// return a map from the representative node of each set to the list of nodes in their sets
///
std::map< platform::Size, utility::vector1< platform::Size > >
DisjointSets::sets() const {

	std::map< platform::Size, utility::vector1< platform::Size > > r2s;

	for ( platform::Size i = 1, ie = n_nodes(); i <= ie; ++i ) {
		// operator[] safe here
		r2s[ ds_find( i ) ].push_back( i );
	}

	return r2s;
}


}
}

