// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// This file contains code derived from the Boost graph library.
// See Rosetta/main/source/external/boost_1_55_0/LICENSE_1_0.txt
// for the Boost library license.

/// @file   utility/graph/BFS_prune.hh
/// @brief  A breadth first search with pruning for boost graphs
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_graph_BFS_prune_HH
#define INCLUDED_utility_graph_BFS_prune_HH

#include <utility/excn/Exceptions.hh>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility/vector1.hh>
#include <platform/types.hh>

namespace utility {
namespace graph {

/// @brief Class to raise to do an immediate stop of a breadth first search.
/// ONLY THROW FROM WITHIN A VISITOR PASSED TO breadth_first_visit_prune/breadth_first_search_prune
class EXCN_Stop_BFS : public utility::excn::Exception {
	using utility::excn::Exception::Exception;
};


/// @brief breadth_first_visit_prune is a slightly
/// modified version of the Boost function breadth_first_visit,
/// allowing the visitor class to prune nodes and edges.
/// See breadth_first_search_prune for details

template <class IncidenceGraph, class Buffer, class BFSVisitor,
class ColorMap>
void breadth_first_visit_prune
(const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor s,
	BFSVisitor vis, ColorMap color, Buffer& Q)
{
	boost::function_requires< boost::IncidenceGraphConcept<IncidenceGraph> >();
	typedef boost::graph_traits<IncidenceGraph> GTraits;
	typedef typename GTraits::vertex_descriptor Vertex;
	//typedef typename GTraits::edge_descriptor Edge;
	boost::function_requires< boost::BFSVisitorConcept<BFSVisitor, IncidenceGraph> >();
	boost::function_requires< boost::ReadWritePropertyMapConcept<ColorMap, Vertex> >();
	typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
	typedef boost::color_traits<ColorValue> Color;
	typename GTraits::out_edge_iterator ei, ei_end;

	try {
		boost::put(color, s, Color::gray());
		if ( vis.discover_vertex(s, g) ) return;
		Q.push(s);
		while ( ! Q.empty() ) {
			Vertex u = Q.top(); Q.pop();
			if ( vis.examine_vertex(u, g) ) continue;
			for ( boost::tie(ei, ei_end) = boost::out_edges(u, g); ei != ei_end; ++ei ) {
				Vertex v = boost::target(*ei, g);
				if ( vis.examine_edge(*ei, g) ) continue;
				ColorValue v_color = get(color, v);
				if ( v_color == Color::white() ) {
					if ( vis.tree_edge(*ei, g) ) continue;
					boost::put(color, v, Color::gray());
					if ( vis.discover_vertex(v, g) ) continue;
					Q.push(v);
				} else {
					if ( vis.non_tree_edge(*ei, g) ) continue;
					if ( v_color == Color::gray() )       vis.gray_target(*ei, g);
					else                                vis.black_target(*ei, g);
				}
			} // end for
			boost::put(color, u, Color::black());          vis.finish_vertex(u, g);
		} // end while
	} catch ( EXCN_Stop_BFS const & e ) {
		; // Do nothing. The exception was just there to cut through the remaining portions.
	}

} // breadth_first_visit

template <class IncidenceGraph, class BFSVisitor,
class ColorMap>
void breadth_first_visit_prune
(const IncidenceGraph &,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor s,
	BFSVisitor vis, ColorMap color)
{
	boost::queue< typename boost::graph_traits<IncidenceGraph>::vertex_descriptor > Q;
	breadth_first_visit_prune( s, vis, color, Q );
}

/// @brief breadth_first_search_prune is a slightly
/// modified versions of the Boost functions breadth_first_search
/// allowing the visitor class to prune nodes and edges.
///
/// Note the calling order is slightly different (and has to be explicit),
/// as I didn't want to recapitulate the boost frontend magic.
///
/// It assumes all of the relevant functions in the visitor class return bools, with the following meanings:
///
///    vis.initialize_vertex(u, g) -- return value ignored
///    vis.discover_vertex(u, g);  -- if true, don't add vertex to queue (but still mark as grey)
///    vis.examine_vertex(u, g);   -- if true, don't examine any out edges for vertex (but still mark black)
///    vis.examine_edge(e, g);     -- if true, ignore edge (don't follow)
///    vis.tree_edge(e, g);        -- if true, ignore discovered vertex (don't mark as grey)
///    vis.non_tree_edge(e, g);    -- if true, don't bother with black/grey function calls
///    vis.gray_target(e, g);      -- return value ignored
///    vis.black_target(e, g);     -- return value ignored
///    vis.finish_vertex(u, g);    -- return value ignored
///
/// Any of the above functions can throw a EXCN_Stop_BFS exception, which will immediately halt the search.

template <class VertexListGraph, class Buffer, class BFSVisitor,
class ColorMap>
void breadth_first_search_prune
(const VertexListGraph& g,
	typename boost::graph_traits<VertexListGraph>::vertex_descriptor s,
	BFSVisitor vis, ColorMap color, Buffer& Q )
{
	// Initialization
	typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
	typedef boost::color_traits<ColorValue> Color;
	typename boost::graph_traits<VertexListGraph>::vertex_iterator i, i_end;

	try {
		for ( boost::tie(i, i_end) = boost::vertices(g); i != i_end; ++i ) {
			vis.initialize_vertex(*i, g);
			boost::put(color, *i, Color::white());
		}
		breadth_first_visit_prune(g, s, vis, color, Q);
	} catch ( EXCN_Stop_BFS const & e ) {
		; // Do nothing. Exception already halted.
	}
}

template <class VertexListGraph, class BFSVisitor, class ColorMap>
void breadth_first_search_prune
(const VertexListGraph& g,
	typename boost::graph_traits<VertexListGraph>::vertex_descriptor s,
	BFSVisitor vis, ColorMap color)
{
	boost::queue< typename boost::graph_traits<VertexListGraph>::vertex_descriptor > Q;
	breadth_first_visit_prune(g, s, vis, color, Q );
}

template <class VertexListGraph, class BFSVisitor>
void breadth_first_search_prune
(const VertexListGraph& g,
	typename boost::graph_traits<VertexListGraph>::vertex_descriptor s,
	BFSVisitor vis)
{
	boost::queue< typename boost::graph_traits<VertexListGraph>::vertex_descriptor > Q;
	breadth_first_visit_prune(g, s, vis, make_two_bit_color_map(num_vertices(g), boost::get(boost::vertex_index,g)), Q );
}

class null_bfs_prune_visitor {
public:
	template <class Vertex, class Graph>
	bool initialize_vertex(Vertex, Graph& ) { return false; }

	template <class Vertex, class Graph>
	bool discover_vertex(Vertex , Graph& ) { return false; }

	template <class Vertex, class Graph>
	bool examine_vertex(Vertex , Graph& ) { return false; }

	template <class Edge, class Graph>
	bool examine_edge(Edge , Graph& ) { return false; }

	template <class Edge, class Graph>
	bool tree_edge(Edge , Graph& ) { return false; }

	template <class Edge, class Graph>
	bool non_tree_edge(Edge , Graph& ) { return false; }

	template <class Edge, class Graph>
	bool gray_target(Edge , Graph& ) { return false; }

	template <class Edge, class Graph>
	bool black_target(Edge , Graph& ) { return false; }

	template <class Vertex, class Graph>
	bool finish_vertex(Vertex , Graph& ) { return false; }

};


template< class Graph>
class HideVertexVisitor: public utility::graph::null_bfs_prune_visitor {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VD;

private:
	// References to avoid copying when the visitor is passed by value.
	VD hidden_vertex_;
	platform::Size & number_valid_vertices_;
	utility::vector1<VD> & connected_vertices_;
public:
	HideVertexVisitor(VD const & hidden, platform::Size & size, utility::vector1<VD> & connected_vertices):
		hidden_vertex_( hidden ),
		number_valid_vertices_( size ),
		connected_vertices_( connected_vertices )
	{}

	platform::Size size(){return number_valid_vertices_;}
	utility::vector1<VD> vertices(){return connected_vertices_;}

	bool discover_vertex(VD vertex, Graph const & /*graph*/) {
		if ( vertex == hidden_vertex_ ) {
			return true;
		} else {
			++number_valid_vertices_;
			connected_vertices_.push_back(vertex);
			return false;
		}
	}

};

} // namespace graph
} // namespace utility

#endif
