// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// @file ResidueGraphTypes.hh
///
/// @brief Graph structure for ResidueType
///
/// @details
/// This is the typedefs for the graph implementation of ResidueType based on boost graphs.
///
/// @author Steven Combs
////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_chemical_ResidueGraphTypes_hh
#define INCLUDED_core_chemical_ResidueGraphTypes_hh

// Unit headers
#include <core/chemical/Atom.hh> // needed full header for ResidueGraph def
#include <core/chemical/Bond.hh>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <boost/graph/adjacency_list.hpp>
#include <utility>

// Package headers

namespace core {
namespace chemical {

/////////////////////////////////////////
/////////// Graph typedefs //////////////

typedef boost::undirected_graph<
	Atom, // struct with properties of a node
	Bond // struct with properties of an edge
	/*,ResidueType*/
	> ResidueGraph;

typedef ResidueGraph::vertex_descriptor VD;
typedef ResidueGraph::edge_descriptor ED;
typedef utility::vector1< VD > VDs;

// A note regarding iterator "const-ness":
// Unlike iterators for the standard library containers, which can be dereferenced to directly obtain container members,
// the boost graph iterators dereference to vertex/edge descriptors, which are functionally equivalent to index values
// (though internal implementation details differ).
// The const-ness of the container doesn't transfer over to the vertex/edge descriptors,
// any more than the const-ness of a utility::vector1 transfers over to the core::Size you use to index it.
// For this reason, you don't need both const and non-const versions of boost::graph iterators,
// and functions returning iterators can be considered 'const' against the class/graph.
// (Actual modification of the graph would require a seperate non-const function call against the graph.)
//
// ref: http://lists.boost.org/Archives/boost/2001/07/14838.php

typedef boost::graph_traits<ResidueGraph>::vertex_iterator VIter;
typedef std::pair<VIter, VIter> VIterPair;

typedef boost::graph_traits<ResidueGraph>::edge_iterator EIter;
typedef std::pair<EIter, EIter> EIterPair;

typedef boost::graph_traits<ResidueGraph>::out_edge_iterator OutEdgeIter;
//typedef boost::graph_traits<ResidueGraph>::in_edge_iterator InEdgeIter; // Out and in edges are the same in an undirected_graph
typedef std::pair<OutEdgeIter, OutEdgeIter> OutEdgeIterPair;

typedef boost::graph_traits<ResidueGraph>::adjacency_iterator AdjacentIter;
typedef std::pair<AdjacentIter, AdjacentIter> AdjacentIterPair;

typedef std::map< std::string, VD > NameVDMap;
typedef std::pair<std::string, VD> NameVDPair;
typedef std::pair<NameVDMap::iterator, bool> NameVDInserted;

////////////////////////////////////////////////
/////////// Convenience Functions //////////////

// These are here because they depend somewhat on the particular implementation choice of ResidueGraph

/// @brief Does a ResidueGraph have a given vertex descriptor?
inline bool
has( ResidueGraph const & graph, VD vd ) {
	VIterPair iters( boost::vertices(graph) );
	return std::find( iters.first, iters.second, vd ) != iters.second;
}
/// @brief Does a ResidueGraph have a given edge descriptor?
inline bool
has( ResidueGraph const & graph, ED ed ) {
	EIterPair iters( boost::edges(graph) );
	return std::find( iters.first, iters.second, ed ) != iters.second;
}

/// @brief When adding and deleting nodes in a graph, sometimes the inner counting of nodes/edges gets outdated.
///Run this to fix the problem.
template <typename Graph>
void regenerate_graph_vertex_index(Graph & graph){
	Size index(0); //counter
	typename boost::graph_traits<Graph>::vertex_iterator itr, end;
	for ( boost::tie(itr,end) = boost::vertices(graph); itr !=end; ++itr ) {
		boost::put(boost::vertex_index, graph, *itr, index);
		++index;
	}
}

///Light weight graph typedefs
///The light weight graph is a graph that holds a pointer to the edge descriptor
///and vertex descriptor or the ResidueGraph. We generate the light weight graph
///so that we can do rapid things like look for rings ina small molecule. Also
//the properties are vectors, so there is fast random lookup
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_name_t, core::chemical::VD>,
	boost::property<boost::edge_name_t, core::chemical::ED> > LightWeightResidueGraph;
typedef boost::graph_traits<LightWeightResidueGraph>::vertex_descriptor lwrg_VD;
typedef boost::graph_traits<LightWeightResidueGraph>::edge_descriptor lwrg_ED;
typedef boost::graph_traits<LightWeightResidueGraph>::vertex_iterator lwrg_vd_iter;
typedef std::pair<lwrg_vd_iter, lwrg_vd_iter> lwrg_vd_pair_iter;
typedef boost::graph_traits<LightWeightResidueGraph>::edge_iterator lwrg_edge_iter;
typedef boost::graph_traits<LightWeightResidueGraph>::out_edge_iterator lwrg_out_edge_iter;
typedef std::pair<lwrg_out_edge_iter, lwrg_out_edge_iter> lwrg_out_edge_iter_pair;

}
}

#endif
