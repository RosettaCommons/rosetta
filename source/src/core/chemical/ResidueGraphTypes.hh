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
/// @brief
/// Graph structure for ResidueType
///
/// @details
/// This is the typedefs and filtered graphs for the graph implementation of ResidueType based on boost graphs.
/// Filtered graphs are graph structures that have been filtered based on a certain criteria. For example, the
/// Acceptor atom graph has been filtered so that every node and edge in the graph is associated with an acceptor
/// atom. The properties of the filtered graphs can be determined by any criteria. Currently, atom types are used
/// as the metric to filter the graphs. This does not have to be the case. Graphs can be filtered based on the
/// atoms, orbitals, etc etc. It is up to your imagination. The unit tests for these show examples of how to use
/// the filtered graphs.
///
/// Each filter graph has an operator that is used to determine if a node should be in the graph. An iterator through
/// each node and edge of the graph is available. Specifically, if you want to iterate through the graph nodes, you would
/// use this method: for(HeavyAtomVIterPair vp = boost::vertices(heavy_atom_graph); vp.first != vp.second; ++vp.first){}
///
/// @author Steven Combs
////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_chemical_ResiduePredicates_hh
#define INCLUDED_core_chemical_ResiduePredicates_hh

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////////// PREDICATES for FILTERED GRAPHS ///////////////////
////////////////////////////////////////////////////////////


/// @default A filtered graph that doesn't contain fake/virtual atoms and fake/virtual bonds.

class RealFilter{
public:
	RealFilter()= default;
	RealFilter(ResidueGraph const & graph):graph_(&graph) {};
	bool operator()(VD const vd) const;
	bool operator()(ED const ed) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
};
typedef boost::filtered_graph<ResidueGraph, RealFilter, RealFilter> RealResidueGraph;
typedef RealResidueGraph::vertex_descriptor RealResidueVD;
typedef RealResidueGraph::edge_descriptor RealResidueED;
typedef boost::graph_traits<RealResidueGraph>::vertex_iterator RealResidueVIter;
typedef boost::graph_traits<RealResidueGraph>::edge_iterator RealResidueEIter;
typedef boost::graph_traits<RealResidueGraph>::out_edge_iterator RealResidueOutEdgeIter;
typedef std::pair<RealResidueOutEdgeIter, RealResidueOutEdgeIter> RealResidueOutEdgeIterPair;
typedef std::pair<RealResidueVIter, RealResidueVIter> RealResidueVIterPair;
typedef boost::graph_traits<RealResidueGraph>::adjacency_iterator RealResidueAdjacentIter;
typedef std::pair<RealResidueAdjacentIter, RealResidueAdjacentIter> RealResidueAdjacentIterPair;

/// @brief The filter responsible for obtaining all heavy atoms.
class HeavyAtomFilter{
public:
	HeavyAtomFilter(){};
	HeavyAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomFilter> HeavyAtomGraph;
typedef HeavyAtomGraph::vertex_descriptor HeavyAtomVD;
typedef HeavyAtomGraph::edge_descriptor HeavyAtomED;
typedef boost::graph_traits<HeavyAtomGraph>::vertex_iterator HeavyAtomVIter;
typedef boost::graph_traits<HeavyAtomGraph>::edge_iterator HeavyAtomEIter;
typedef boost::graph_traits<HeavyAtomGraph>::out_edge_iterator HeavyAtomOutEdgeIter;
typedef std::pair<HeavyAtomOutEdgeIter, HeavyAtomOutEdgeIter> HeavyAtomOutEdgeIterPair;
typedef std::pair<HeavyAtomVIter, HeavyAtomVIter> HeavyAtomVIterPair;


////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for obtaining all acceptor atoms.
class AcceptorAtomFilter{
public:
	AcceptorAtomFilter(){};
	AcceptorAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, AcceptorAtomFilter> AcceptorAtomGraph;
typedef boost::graph_traits<AcceptorAtomGraph>::vertex_iterator AcceptorAtomVIter;
typedef boost::graph_traits<AcceptorAtomGraph>::edge_iterator AcceptorAtomEIter;
typedef std::pair<AcceptorAtomVIter, AcceptorAtomVIter> AcceptorAtomVIterPair;

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for obtaining all heavy atoms with polar hydrogens attached to them.
class HeavyAtomWithPolarHydrogensFilter{
public:
	HeavyAtomWithPolarHydrogensFilter(){};
	HeavyAtomWithPolarHydrogensFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomWithPolarHydrogensFilter> HeavyAtomWithPolarHydrogensGraph;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::vertex_iterator HeavyAtomWithPolarHydrogensVIter;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::edge_iterator HeavyAtomWithPolarHydrogensEIter;
typedef std::pair<HeavyAtomWithPolarHydrogensVIter, HeavyAtomWithPolarHydrogensVIter> HeavyAtomWithPolarHydrogensVIterPair;

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for finding heavy atoms with hydrogens.
class HeavyAtomWithHydrogensFilter{
public:
	HeavyAtomWithHydrogensFilter(){};
	HeavyAtomWithHydrogensFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomWithHydrogensFilter> HeavyAtomWithHydrogensGraph;
typedef boost::graph_traits<HeavyAtomWithHydrogensGraph>::vertex_iterator HeavyAtomWithHydrogensVIter;
typedef boost::graph_traits<HeavyAtomWithHydrogensGraph>::edge_iterator HeavyAtomWithHydrogensEIter;
typedef std::pair<HeavyAtomWithHydrogensVIter, HeavyAtomWithHydrogensVIter> HeavyAtomWithHydrogensVIterPair;

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


/// @brief The filter responsible for all hydrogens.
class HydrogenAtomFilter{
public:
	HydrogenAtomFilter(){};
	HydrogenAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HydrogenAtomFilter> HydrogenAtomGraph;
typedef HydrogenAtomGraph::vertex_descriptor HydrogenAtomVD;
typedef HydrogenAtomGraph::edge_descriptor HydrogenAtomED;
typedef boost::graph_traits<HydrogenAtomGraph>::vertex_iterator HydrogenAtomVIter;
typedef boost::graph_traits<HydrogenAtomGraph>::edge_iterator HHydrogenAtomEIter;
typedef boost::graph_traits<HydrogenAtomGraph>::out_edge_iterator HydrogenAtomOutEdgeIter;
typedef std::pair<HydrogenAtomOutEdgeIter, HydrogenAtomOutEdgeIter> HydrogenAtomOutEdgeIterPair;
typedef std::pair<HydrogenAtomVIter, HydrogenAtomVIter> HydrogenAtomVIterPair;


////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for all aromatic atoms.
class AromaticAtomFilter{
public:
	AromaticAtomFilter(){};
	AromaticAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, AromaticAtomFilter> AromaticAtomGraph;
typedef boost::graph_traits<AromaticAtomGraph>::vertex_iterator AromaticAtomVIter;
typedef boost::graph_traits<AromaticAtomGraph>::edge_iterator AromaticAtomEIter;
typedef std::pair<AromaticAtomVIter, AromaticAtomVIter> AromaticAtomVIterPair;


////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for all polar hydrogens.
class PolarHydrogenFilter{
public:
	PolarHydrogenFilter(){};
	PolarHydrogenFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, PolarHydrogenFilter> PolarHydrogenGraph;
typedef boost::graph_traits<PolarHydrogenGraph>::vertex_iterator PolarHydrogenVIter;
typedef boost::graph_traits<PolarHydrogenGraph>::edge_iterator PolarHydrogenEIter;
typedef std::pair<PolarHydrogenVIter, PolarHydrogenVIter> PolarHydrogenVIterPair;

/// @brief The filter responsible for all apolar hydrogens.
class APolarHydrogenFilter{
public:
	APolarHydrogenFilter(){};
	APolarHydrogenFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, APolarHydrogenFilter> APolarHydrogenGraph;
typedef boost::graph_traits<APolarHydrogenGraph>::vertex_iterator APolarHydrogenVIter;
typedef boost::graph_traits<APolarHydrogenGraph>::edge_iterator APolarHydrogenEIter;
typedef std::pair<APolarHydrogenVIter, APolarHydrogenVIter> APolarHydrogenVIterPair;


template< class Graph1, class Graph2 >
class CopyVertex {
public:
	CopyVertex(const Graph1& g1, Graph2& g2):
		g1_( g1 ),
		g2_(g2 )
	{ }

	void operator()(const typename Graph1::vertex_descriptor& v1, typename Graph2::vertex_descriptor& v2) const {
		g2_[v2] = g1_[v1];
	}

private:
	Graph1 const & g1_;
	Graph2 & g2_;

};


template< class Graph1, class Graph2 >
class CopyEdge {
public:
	CopyEdge(const Graph1& g1, Graph2& g2):
		g1_( g1 ),
		g2_(g2 )
	{ }

	void operator()(const typename Graph1::edge_descriptor& e1, typename Graph2::edge_descriptor& e2) const {
		g2_[e2] = g1_[e1];
	}
private:

	Graph1 const & g1_;
	Graph2 & g2_;

};

}
}

#endif
