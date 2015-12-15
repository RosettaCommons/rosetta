 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/UpperEdgeGraph.srlz.hh
/// @brief  templated graph for fast edge additions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef    SERIALIZATION

#ifndef INCLUDED_core_graph_UpperEdgeGraph_SRLZ_HH
#define INCLUDED_core_graph_UpperEdgeGraph_SRLZ_HH

// Unit Headers
#include <core/graph/UpperEdgeGraph.hh>

// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>


namespace core {
namespace graph {

template < class Archive, class V, class E >
void save_to_archive( Archive & arc, typename core::graph::UpperEdgeGraph< V, E > const & graph ) {
	arc( graph.num_vertices() );
	for ( Size ii = 1; ii <= graph.num_vertices(); ++ii ) {
		arc( graph.get_vertex( ii ).data() );
	}
	for ( Size ii = 1; ii <= graph.num_vertices(); ++ii ) {
		UEVertex< V, E > const & iivert( graph.get_vertex( ii ) );
		arc( iivert.num_upper_neighbors() );
		for ( typename UEVertex< V, E >::UpperEdgeListConstIter
				iter = iivert.const_upper_edge_list_begin(), iter_end = iivert.const_upper_edge_list_end();
				iter != iter_end; ++iter ) {
			if ( iter->deleted() ) { continue; }
			arc( iter->upper_vertex() );
			arc( iter->data() );
		}
	}
}

template < class Archive, class V, class E >
void load_from_archive( Archive & arc, typename core::graph::UpperEdgeGraph< V, E > & graph ) {
	platform::Size nvertices( 0 ); arc( nvertices );
	graph.set_num_vertices( nvertices );
	for ( Size ii = 1; ii <= graph.num_vertices(); ++ii ) {
		arc( graph.get_vertex( ii ).data() );
	}
	for ( Size ii = 1; ii <= graph.num_vertices(); ++ii ) {
		platform::Size ii_n_upper_neighbors( 0 ); arc( ii_n_upper_neighbors );
		for ( Size jj = 1; jj <= ii_n_upper_neighbors; ++jj ) {
			platform::Size jj_upper_neighbor( 0 ); arc( jj_upper_neighbor );
			E edge_data;
			arc( edge_data );
			graph.add_edge( ii, jj_upper_neighbor, edge_data );
		}
	}
}

}
}

#endif

#endif // SERIALIZATION
