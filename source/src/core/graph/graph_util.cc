// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/graph_util.hh
/// @brief  generic graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/graph/graph_util.hh>

// Package Headers
#include <core/graph/Graph.hh>

//#include <iostream>

namespace core {
namespace graph {

bool
all_visited( utility::vector1< bool > const & visited );

/// @details DFS search to identify connected components
/// O(V) memory, O(V+E) time.
utility::vector1< std::pair< platform::Size, platform::Size > >
find_connected_components( Graph const & g )
{
	using namespace utility;
	vector1< bool > reached( g.num_nodes(), false );

	/// CC description
	vector1< platform::Size > representative;
	vector1< platform::Size > cc_nelements;
	representative.reserve( g.num_nodes() ); // O(N) -- at most N cc's.
	cc_nelements.reserve( g.num_nodes() );   // O(N) -- at most N cc's.

	/// DFS data
	vector1< platform::Size > to_explore( g.num_nodes(), 0 ); // insert at most N
	platform::Size n_to_explore = 0;

	for ( platform::Size ii = 1; ii <= (platform::Size) g.num_nodes(); ++ii ) {
		if ( reached[ ii ] ) continue;

		// new connected component
		// run a dfs from this node
		representative.push_back( ii );
		cc_nelements.push_back( 1 );
		n_to_explore = 0;
		to_explore[ ++n_to_explore ] = ii;

		while ( n_to_explore != 0 ) {
			platform::Size const exploring = to_explore[ n_to_explore ];
			to_explore[ n_to_explore ] = 0; // for sanity -- unneccesary
			--n_to_explore;

			Node const * node_exploring = g.get_node( exploring );
			for ( EdgeListConstIterator
					eiter = node_exploring->const_edge_list_begin(),
					eiter_end = node_exploring->const_edge_list_end();
					eiter != eiter_end; ++eiter ) {
				platform::Size neighbor = (*eiter)->get_other_ind( exploring );
				if ( ! reached[ neighbor ] ) {
					to_explore[ ++n_to_explore ] = neighbor; // dfs
					++cc_nelements[ cc_nelements.size() ];
					reached[ neighbor ] = true;
				}
			}
		}
	}

debug_assert( all_visited( reached ) );

	// Prepare output descriptions
	vector1< std::pair< platform::Size, platform::Size > > cc_descriptions( representative.size() );
	for ( platform::Size ii = 1; ii <= representative.size(); ++ii ) {
		cc_descriptions[ ii ] = std::make_pair( representative[ ii ], cc_nelements[ ii ] );
	}
	return cc_descriptions;
}

bool
all_visited( utility::vector1< bool > const & visited )
{
	for ( platform::Size ii = 1; ii <= visited.size(); ++ii ) {
		if ( ! visited[ ii ] ) return false;
	}
	return true;
}

utility::vector1< std::pair< platform::Size, platform::Size > >
find_connected_components( Graph const & g );

void
delete_all_intragroup_edges(
	Graph & g,
	utility::vector1< platform::Size > const & node_groups
)
{
debug_assert( node_groups.size() == g.num_nodes() );

	for( Graph::EdgeListIter edge_it = g.edge_list_begin();	 edge_it != g.edge_list_end(); ){

		Graph::EdgeListIter next_edge = edge_it;
		++next_edge;

		if( node_groups[ (*edge_it)->get_first_node_ind() ] == node_groups[  (*edge_it)->get_second_node_ind() ] ){
			//std::cout << "GRAPH_EDGE_DELETE resi " << (*edge_it)->get_first_node_ind() << " to " << (*edge_it)->get_second_node_ind() << std::endl;
			g.delete_edge( *edge_it );

		}

		edge_it = next_edge;

	} //loop over edges
} //delete_all_group1_edges_except_to_group2

}
}
