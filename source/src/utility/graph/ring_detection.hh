// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/graph/ring_detection.hh
/// @brief  Algorithms for working with rings in boost graphs
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_graph_ring_detection_HH
#define INCLUDED_utility_graph_ring_detection_HH

#include <utility/graph/BFS_prune.hh>
#include <utility/vector0.hh>
#include <utility/exit.hh>

#include <platform/types.hh>

#include <boost/graph/undirected_dfs.hpp>

#include <map>
#include <set>

namespace utility {
namespace graph {

/// @brief A class to implement the behavior of the smallest ring size finding algorithm, accessible through the smallest_ring_size() function below.
/// @details Based on BCL's smallest ring size detection algorithm.

template< class Graph, class DistanceMap, class LabelMap >
class RingSizeVisitor: public utility::graph::null_bfs_prune_visitor {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VD;

private:
	// References to avoid copying when the visitor is passed by value.
	DistanceMap & distances_;
	LabelMap & labels_;
	// Reference to allow output when the visitor is passed by value
	platform::Size & size_;
	platform::Size stop_level_;
public:
	RingSizeVisitor(VD const & source, Graph const & graph, DistanceMap & distances, LabelMap & labels, platform::Size & size, platform::Size const & max_size = 2*999999):
		distances_( distances ),
		labels_( labels ),
		size_(size), // Tie for output.
		stop_level_( max_size/2 + 1 ) // Integer truncation division intended
	{
		size_ = 999999;
		boost::put(distances_,source,0);
		boost::put(labels_,source,0);
		typename boost::graph_traits<Graph>::adjacency_iterator iter, end;
		platform::Size index=1;
		for( boost::tie(iter,end) = boost::adjacent_vertices(source,graph); iter != end; ++iter, ++index) {
			boost::put(labels_,*iter,index);
			boost::put(distances_,*iter,1);
		}
	}

	bool tree_edge(Edge const & e, Graph const & g) {
		//First time we see the target of the edge - check validity, mark distance and label
		// We can assume the parent node is initialized on the data map.
		VD const & parent( boost::source(e,g) );
		VD const & child( boost::target(e,g) );
		platform::Size distance( boost::get(distances_,parent)+1 );
		if( distance == 1 ) { return false; } // We've already processed the node in the constructor.
		if( distance >= stop_level_ ) { return true; }

		//Set Distance and label
		boost::put(distances_, child, distance);
		boost::put(labels_, child, boost::get(labels_,parent));
		return false;
	}

	bool gray_target(Edge const & e, Graph const & g) {
		// A grey target implies a ring closure.
		// We can be assured that both nodes have been initialized on the data maps.
		VD const & left( boost::source(e,g) );
		VD const & right( boost::target(e,g) );
		if( boost::get(labels_,left) ==  boost::get(labels_,right) ) {
			// We're closing a ring that can bypass the start vertex. Ignore.
			return false;
		}
		platform::Size const & dleft( boost::get(distances_,left) );
		platform::Size const & dright( boost::get(distances_,right) );
		platform::Size ringsize( dleft + dright + 1 ); // +1 for the start vertex
		// Note that the first ring encountered may not be the smallest ring.
		// E.g. if we're working on the second level, we may close a six member ring first,
		// but then close a five member one later while finishing up the second level.
		// Finish this level, but don't bother working on the next level.
		if( ringsize < size_ ) {
			size_ = ringsize;
			stop_level_ = dleft + 1; // Don't follow any nodes to next level
		}
		return false;
	}
};

/// @brief A boost graph-based function to find the size of the smallest ring for a given vertex.
/// Will return the maximum ring size, or 999999 for no ring. If there is a maximum ring size,
/// That can be set to limit the amount of search needed.
template < class Graph >
platform::Size
smallest_ring_size( typename boost::graph_traits< Graph>::vertex_descriptor const & vd, Graph const & graph, platform::Size const & max_ring_size = 2*999999 ) {

	typedef typename boost::graph_traits< Graph>::vertex_descriptor VD;

/*	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIDMap;
	//vector0 as the VertexIDMap is zero-based.
	typedef utility::vector0< platform::Size > DistanceMapStore;
	typedef boost::iterator_property_map<DistanceMapStore::iterator, VertexIDMap,
			std::iterator_traits<DistanceMapStore::iterator>::value_type,
			std::iterator_traits<DistanceMapStore::iterator>::reference
			> DistanceMap;
	typedef utility::vector0< platform::Size > LabelMapStore;
	typedef boost::iterator_property_map<LabelMapStore::iterator, VertexIDMap,
			std::iterator_traits<LabelMapStore::iterator>::value_type,
			std::iterator_traits<LabelMapStore::iterator>::reference
			> LabelMap;

	VertexIDMap vertex_id( boost::get(boost::vertex_index, graph) );
	DistanceMapStore	distancestore( boost::num_vertices(graph), 999999 );
	DistanceMap	distances( boost::make_iterator_property_map(distancestore.begin(), vertex_id) );
	LabelMapStore	labelstore( boost::num_vertices(graph), 999999 );
	LabelMap	labels( boost::make_iterator_property_map(labelstore.begin(), vertex_id) );
*/

	typedef typename std::map<VD,platform::Size> DataStore;
	typedef typename boost::associative_property_map< DataStore > DataStoreMap;

	DataStore distancestore;
	DataStoreMap distances(distancestore);
	DataStore labelstore;
	DataStoreMap labels(labelstore);

	platform::Size smallest_size = 999999;
	RingSizeVisitor< Graph, DataStoreMap, DataStoreMap > vis(vd, graph, distances, labels, smallest_size, max_ring_size);
	utility::graph::breadth_first_search_prune(graph, vd, vis);

	return smallest_size;
}


template< class Graph, class EdgeMap, class PathMap >
class RingEdgeAnnotationVisitor: public boost::default_dfs_visitor {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VD;

private:
	// References to avoid copying when the visitor is passed by value.
	// EdgeMap is a map of VD:(map of VD:bools), where EdgeMap[vd1][vd2] is true if the vd1-vd2 edge is in a cycle.
	EdgeMap & edgemap_;
	// PathMap is a map of VertexDescriptors:set<edge_descriptors>, representing the initial path to the vertex
	PathMap & pathmap_;
public:
	RingEdgeAnnotationVisitor(EdgeMap & edgemap, PathMap & pathmap):
		edgemap_( edgemap ),
		pathmap_( pathmap )
	{}

	void examine_edge(Edge, const Graph&) {
		//edgemap_[u] = false; // Initialize edge as non-cycle.
		// // Maps default to zero initialized if not found.
	}
	void tree_edge(Edge u, const Graph& g) {
		//std::cout << "Tree: " << boost::source(u,g) << "---" << boost::target(u,g) << std::endl;
		typename PathMap::mapped_type path( pathmap_[boost::source(u,g)] ); // make a copy.
		path.insert( u );
		pathmap_[ boost::target(u,g) ] = path;
	}
	void back_edge(Edge u, const Graph& g) {
		//std::cout << "Back: " << boost::source(u,g) << "---" << boost::target(u,g) << std::endl;
		typename PathMap::mapped_type const & path( pathmap_[boost::source(u,g)] );
		typename PathMap::mapped_type const & path2( pathmap_[boost::target(u,g)] );
		// We've found a cycle
		edgemap_[boost::source(u,g)][boost::target(u,g)] = true; //The back edge is always part of the cycle.
		edgemap_[boost::target(u,g)][boost::source(u,g)] = true;
		typename PathMap::mapped_type::const_iterator iter, iter_end;
		for( iter = path.begin(), iter_end = path.end(); iter != iter_end; ++iter ) {
			if( path2.count( *iter ) == 0 ) {
				// If we have an edge in the longer path, but it isn't on the "stem"
				// to the node we're going back to, it's part of a cycle.
				//std::cout << "Marking Ring " << boost::source(*iter,g) << "---" << boost::target(*iter,g) << std::endl;
				edgemap_[boost::source(*iter,g)][boost::target(*iter,g)] = true;
				edgemap_[boost::target(*iter,g)][boost::source(*iter,g)] = true;
			}
		}
	}
	void forward_or_cross_edge(Edge, const Graph&) {
		// Have to use generic tracer - can't initialize static tracer in header file.
		//std::cout << "Cross: " << boost::source(u,g) << "---" << boost::target(u,g) << std::endl;
		utility_exit_with_message( "Found forward or cross edge when detecting cycles - should not happen with an undirected graph." );
	}
};

template < class Graph >
std::map< typename boost::graph_traits<Graph>::vertex_descriptor, std::map< typename boost::graph_traits<Graph>::vertex_descriptor, bool > >
annotate_ring_edges( Graph & graph) {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VD;
	typedef typename std::map< typename boost::graph_traits<Graph>::vertex_descriptor, std::map< typename boost::graph_traits<Graph>::vertex_descriptor, bool > > EdgeMap;
	typedef typename std::map< VD, std::set< Edge > >  PathMap;
	typedef typename std::map< VD, boost::default_color_type> VertexColorMap;
	typedef typename std::map< Edge, boost::default_color_type> EdgeColorMap;

	EdgeMap edgemap;
	PathMap pathmap;

	VertexColorMap vertexcolormap;
	EdgeColorMap edgecolormap;
	boost::associative_property_map<VertexColorMap> vertexcolorproperty(vertexcolormap);
	boost::associative_property_map<EdgeColorMap> edgecolorproperty(edgecolormap);

	RingEdgeAnnotationVisitor<Graph, EdgeMap, PathMap > vis( edgemap, pathmap );
	//Use DFS here so that we only have to keep track of one active path
	//Use undirected_dfs() because the regular DFS has a back edge for each forward edge.
	boost::undirected_dfs(graph, vis, vertexcolorproperty, edgecolorproperty);

	return edgemap;
}

} // namespace graph
} // namespace utility

#endif
