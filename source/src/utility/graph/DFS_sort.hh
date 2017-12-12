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

/// @file   utility/graph/DFS_sort.hh
/// @brief  A depth first search with sorting for edge visitation.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_graph_DFS_sort_HH
#define INCLUDED_utility_graph_DFS_sort_HH

#include <boost/graph/depth_first_search.hpp>

namespace utility {
namespace graph {


/// @brief depth_first_visit_sort_impl is a slightly
/// modified version of the recursive version of the
/// Boost function depth_first_visit_impl,
/// allowing the visitor class to prune nodes and edges.
/// See depth_first_search_sort for details

template <class IncidenceGraph, class DFSVisitor, class ColorMap,
	class SortFunc>
void depth_first_visit_sort_impl
(const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor u,
	DFSVisitor& vis,  // pass-by-reference here, important!
	ColorMap color, SortFunc const & func)
{
	using namespace boost;
	typedef typename graph_traits<IncidenceGraph>::vertex_descriptor Vertex;
	typedef typename graph_traits<IncidenceGraph>::edge_descriptor Edge;
	typedef std::vector< Edge > Edges;
	typedef typename property_traits<ColorMap>::value_type ColorValue;
	typedef color_traits<ColorValue> Color;

	put(color, u, Color::gray());          vis.discover_vertex(u, g);

	// typedef typename unwrap_reference<SortFunc>::type SF;
	// // Variable is needed to workaround a borland bug.
	// SF& fn = static_cast<SF&>(func);
	Edges edges;
	{ // Anonymous namespace to limit ei/ei_end scope.
		typename graph_traits<IncidenceGraph>::out_edge_iterator ei, ei_end;
		for ( boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
			edges.push_back( *ei );
		}
	}
	std::sort(edges.begin(),edges.end(),func); // Sort edges
	for ( typename Edges::iterator ei( edges.begin() ), ei_end( edges.end() ); ei != ei_end; ++ei ) {
		Vertex v = target(*ei, g);           vis.examine_edge(*ei, g);
		ColorValue v_color = get(color, v);
		if ( v_color == Color::white() ) {
			vis.tree_edge(*ei, g);
			depth_first_visit_sort_impl(g, v, vis, color, func);
		} else if ( v_color == Color::gray() ) {
			vis.back_edge(*ei, g);
		} else {
			vis.forward_or_cross_edge(*ei, g);
		}
	}
	put(color, u, Color::black());         vis.finish_vertex(u, g);
}

/// @brief A sorted depth first search.
/// The parameter func should be a callable object which takes two graph edge descriptors,
/// and returns true if the first should go before the second in ordering.
///
/// This version will start at start_vertex, but will restart arbitrarily for any disconnected graph portions.

template <class VertexListGraph, class DFSVisitor, class ColorMap, class SortFunc>
void
depth_first_search_sort(const VertexListGraph& g, DFSVisitor vis, ColorMap color,
	typename boost::graph_traits<VertexListGraph>::vertex_descriptor start_vertex,
	SortFunc const & func)
{
	using namespace boost;
	typedef typename graph_traits<VertexListGraph>::vertex_descriptor Vertex;
	function_requires<DFSVisitorConcept<DFSVisitor, VertexListGraph> >();
	typedef typename property_traits<ColorMap>::value_type ColorValue;
	typedef color_traits<ColorValue> Color;

	typename graph_traits<VertexListGraph>::vertex_iterator ui, ui_end;
	for ( boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui ) {
		Vertex u = implicit_cast<Vertex>(*ui);
		put(color, u, Color::white()); vis.initialize_vertex(u, g);
	}

	if ( start_vertex != implicit_cast<Vertex>(*vertices(g).first) ) { vis.start_vertex(start_vertex, g);
		depth_first_visit_sort_impl(g, start_vertex, vis, color, func);
	}

	for ( boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui ) {
		Vertex u = implicit_cast<Vertex>(*ui);
		ColorValue u_color = get(color, u);
		if ( u_color == Color::white() ) {       vis.start_vertex(u, g);
			depth_first_visit_sort_impl(g, u, vis, color, func);
		}
	}
}

/// @brief A sorted depth first search.
/// The parameter func should be a callable object which takes two graph edge descriptors,
/// and returns true if the first should go before the second in ordering.
///
/// This version will start at an arbitrary vertex, and restart for disconnected portions of the graph

template <class VertexListGraph, class DFSVisitor, class ColorMap, class SortFunc>
void
depth_first_search_sort(const VertexListGraph& g, DFSVisitor vis, ColorMap color, SortFunc const & func)
{
	using namespace boost;
	typedef typename boost::graph_traits<VertexListGraph>::vertex_iterator vi;
	std::pair<vi, vi> verts = vertices(g);
	if ( verts.first == verts.second ) {
		return;
	}

	depth_first_search_sort(g, vis, color, *verts.first, func);
}

/// @brief Named Parameter Variant
template <class VertexListGraph, class SortFunc, class P, class T, class R>
void
depth_first_search_sort(const VertexListGraph& g, SortFunc const & func,
	const boost::bgl_named_params<P, T, R>& params)
{
	using namespace boost;
	typedef typename boost::graph_traits<VertexListGraph>::vertex_iterator vi;
	std::pair<vi, vi> verts = vertices(g);
	if ( verts.first == verts.second ) {
		return;
	}
	using namespace boost::graph::keywords;
	typedef bgl_named_params<P, T, R> params_type;
	BOOST_GRAPH_DECLARE_CONVERTED_PARAMETERS(params_type, params)
		depth_first_search_sort
		(g,
		arg_pack[_visitor | make_dfs_visitor(null_visitor())],
		boost::detail::make_color_map_from_arg_pack(g, arg_pack),
		arg_pack[_root_vertex | *vertices(g).first],
		func
	);
}

template <class IncidenceGraph, class DFSVisitor, class ColorMap, class SortFunc>
void depth_first_visit
(const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor u,
	DFSVisitor vis, ColorMap color, SortFunc const & func)
{
	vis.start_vertex(u, g);
	depth_first_visit_sort_impl(g, u, vis, color, func);
}

template <class IncidenceGraph, class DFSVisitor, class ColorMap, class SortFunc>
void depth_first_visit
(const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor u,
	DFSVisitor vis, ColorMap color, SortFunc func = SortFunc())
{
	vis.start_vertex(u, g);
	depth_first_visit_sort_impl(g, u, vis, color, func);
}


} // namespace graph
} // namespace utility

#endif
