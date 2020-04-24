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
/// @file ResidueSubGraphTypes.hh
///
/// @brief Filtered Graphs for ResidueType
///
/// @details
/// This is filtered graphs for the graph implementation of ResidueType based on boost graphs.
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
/// @author Rocco Moretti (rmorettiase@gmail.com)
////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_chemical_ResidueSubGraphTypes_hh
#define INCLUDED_core_chemical_ResidueSubGraphTypes_hh

// Unit headers
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/MutableResidueType.hh>
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

/// @brief It's a frustrating truth, but num_vertices doesn't give the number of filtered vertices
/// for filtered graphs, but instead gives the number of underlying vertices ... this gets the true
/// number of filtered vertices.
template< class Graph >
core::Size filtered_num_vertices( Graph const & graph ) {
	typedef typename Graph::vertex_iterator Viter;
	typedef typename std::pair< Viter, Viter > Viterpair;
	core::Size node_count(0);
	for ( Viterpair vip( boost::vertices(graph) );
			vip.first != vip.second;
			++vip.first, ++node_count ) {
		//std::cout << node_count << " " << graph[ *vip.first ].name() << std::endl;
	}
	return node_count;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////////// PREDICATES for FILTERED GRAPHS ///////////////////
////////////////////////////////////////////////////////////

/// @default A filtered graph that doesn't contain fake/virtual atoms and fake/virtual bonds.

class RealFilter{
public:
	RealFilter(): graph_(nullptr) {}; // This can't be private, because various (unused) default objects need it.
	RealFilter(ResidueGraph const & graph): graph_( &graph ) {};
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

/// @brief Convenience constructor
inline
RealResidueGraph
make_real_residue_graph(MutableResidueType const & restype) {
	RealFilter filter(restype.graph());
	RealResidueGraph fg(restype.graph(), filter, filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for obtaining all heavy atoms.
class HeavyAtomFilter{
public:
	HeavyAtomFilter(){};
	HeavyAtomFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomFilter> HeavyAtomGraph;
typedef HeavyAtomGraph::vertex_descriptor HeavyAtomVD;
typedef HeavyAtomGraph::edge_descriptor HeavyAtomED;
typedef boost::graph_traits<HeavyAtomGraph>::vertex_iterator HeavyAtomVIter;
typedef boost::graph_traits<HeavyAtomGraph>::edge_iterator HeavyAtomEIter;
typedef boost::graph_traits<HeavyAtomGraph>::out_edge_iterator HeavyAtomOutEdgeIter;
typedef std::pair<HeavyAtomOutEdgeIter, HeavyAtomOutEdgeIter> HeavyAtomOutEdgeIterPair;
typedef std::pair<HeavyAtomVIter, HeavyAtomVIter> HeavyAtomVIterPair;

/// @brief Convenience constructor
inline
HeavyAtomGraph
make_heavy_atom_graph(MutableResidueType const & restype) {
	HeavyAtomFilter filter(restype.graph(), restype.atom_type_set_ptr());
	HeavyAtomGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for obtaining all acceptor atoms.
class AcceptorAtomFilter{
public:
	AcceptorAtomFilter(){};
	AcceptorAtomFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, AcceptorAtomFilter> AcceptorAtomGraph;
typedef boost::graph_traits<AcceptorAtomGraph>::vertex_iterator AcceptorAtomVIter;
typedef boost::graph_traits<AcceptorAtomGraph>::edge_iterator AcceptorAtomEIter;
typedef std::pair<AcceptorAtomVIter, AcceptorAtomVIter> AcceptorAtomVIterPair;

/// @brief Convenience constructor
inline
AcceptorAtomGraph
make_acceptor_atom_graph(MutableResidueType const & restype) {
	AcceptorAtomFilter filter(restype.graph(), restype.atom_type_set_ptr());
	AcceptorAtomGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for obtaining all heavy atoms with polar hydrogens attached to them.
class HeavyAtomWithPolarHydrogensFilter{
public:
	HeavyAtomWithPolarHydrogensFilter(){};
	HeavyAtomWithPolarHydrogensFilter(ResidueGraph const & graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomWithPolarHydrogensFilter> HeavyAtomWithPolarHydrogensGraph;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::vertex_iterator HeavyAtomWithPolarHydrogensVIter;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::edge_iterator HeavyAtomWithPolarHydrogensEIter;
typedef std::pair<HeavyAtomWithPolarHydrogensVIter, HeavyAtomWithPolarHydrogensVIter> HeavyAtomWithPolarHydrogensVIterPair;

/// @brief Convenience constructor
inline
HeavyAtomWithPolarHydrogensGraph
make_heavy_atom_with_polar_hydrogens_graph(MutableResidueType const & restype) {
	HeavyAtomWithPolarHydrogensFilter filter(restype.graph(), restype.atom_type_set_ptr());
	HeavyAtomWithPolarHydrogensGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for finding heavy atoms with hydrogens.
class HeavyAtomWithHydrogensFilter{
public:
	HeavyAtomWithHydrogensFilter(){};
	HeavyAtomWithHydrogensFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomWithHydrogensFilter> HeavyAtomWithHydrogensGraph;
typedef boost::graph_traits<HeavyAtomWithHydrogensGraph>::vertex_iterator HeavyAtomWithHydrogensVIter;
typedef boost::graph_traits<HeavyAtomWithHydrogensGraph>::edge_iterator HeavyAtomWithHydrogensEIter;
typedef std::pair<HeavyAtomWithHydrogensVIter, HeavyAtomWithHydrogensVIter> HeavyAtomWithHydrogensVIterPair;

/// @brief Convenience constructor
inline
HeavyAtomWithHydrogensGraph
make_heavy_atom_with_hydrogens_graph(MutableResidueType const & restype) {
	HeavyAtomWithHydrogensFilter filter(restype.graph(), restype.atom_type_set_ptr());
	HeavyAtomWithHydrogensGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


/// @brief The filter responsible for all hydrogens.
class HydrogenAtomFilter{
public:
	HydrogenAtomFilter(){};
	HydrogenAtomFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HydrogenAtomFilter> HydrogenAtomGraph;
typedef HydrogenAtomGraph::vertex_descriptor HydrogenAtomVD;
typedef HydrogenAtomGraph::edge_descriptor HydrogenAtomED;
typedef boost::graph_traits<HydrogenAtomGraph>::vertex_iterator HydrogenAtomVIter;
typedef boost::graph_traits<HydrogenAtomGraph>::edge_iterator HHydrogenAtomEIter;
typedef boost::graph_traits<HydrogenAtomGraph>::out_edge_iterator HydrogenAtomOutEdgeIter;
typedef std::pair<HydrogenAtomOutEdgeIter, HydrogenAtomOutEdgeIter> HydrogenAtomOutEdgeIterPair;
typedef std::pair<HydrogenAtomVIter, HydrogenAtomVIter> HydrogenAtomVIterPair;

/// @brief Convenience constructor
inline
HydrogenAtomGraph
make_hydrogen_atom_graph(MutableResidueType const & restype) {
	HydrogenAtomFilter filter(restype.graph(), restype.atom_type_set_ptr());
	HydrogenAtomGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for all aromatic atoms.
class AromaticAtomFilter{
public:
	AromaticAtomFilter(){};
	AromaticAtomFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, AromaticAtomFilter> AromaticAtomGraph;
typedef boost::graph_traits<AromaticAtomGraph>::vertex_iterator AromaticAtomVIter;
typedef boost::graph_traits<AromaticAtomGraph>::edge_iterator AromaticAtomEIter;
typedef std::pair<AromaticAtomVIter, AromaticAtomVIter> AromaticAtomVIterPair;

/// @brief Convenience constructor
inline
AromaticAtomGraph
make_aromatic_atom_graph(MutableResidueType const & restype) {
	AromaticAtomFilter filter(restype.graph(), restype.atom_type_set_ptr());
	AromaticAtomGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// @brief The filter responsible for all polar hydrogens.
class PolarHydrogenFilter{
public:
	PolarHydrogenFilter(){};
	PolarHydrogenFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, PolarHydrogenFilter> PolarHydrogenGraph;
typedef boost::graph_traits<PolarHydrogenGraph>::vertex_iterator PolarHydrogenVIter;
typedef boost::graph_traits<PolarHydrogenGraph>::edge_iterator PolarHydrogenEIter;
typedef std::pair<PolarHydrogenVIter, PolarHydrogenVIter> PolarHydrogenVIterPair;

/// @brief Convenience constructor
inline
PolarHydrogenGraph
make_polar_hydrogen_graph(MutableResidueType const & restype) {
	PolarHydrogenFilter filter(restype.graph(), restype.atom_type_set_ptr());
	PolarHydrogenGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
/// @brief The filter responsible for all apolar hydrogens.
class APolarHydrogenFilter{
public:
	APolarHydrogenFilter(){};
	APolarHydrogenFilter(ResidueGraph const & graph, AtomTypeSetCOP atom_types):graph_(&graph),atom_types_(std::move(atom_types)){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph const * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCOP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, APolarHydrogenFilter> APolarHydrogenGraph;
typedef boost::graph_traits<APolarHydrogenGraph>::vertex_iterator APolarHydrogenVIter;
typedef boost::graph_traits<APolarHydrogenGraph>::edge_iterator APolarHydrogenEIter;
typedef std::pair<APolarHydrogenVIter, APolarHydrogenVIter> APolarHydrogenVIterPair;

/// @brief Convenience constructor
inline
APolarHydrogenGraph
make_apolar_hydrogen_graph(MutableResidueType const & restype) {
	APolarHydrogenFilter filter(restype.graph(), restype.atom_type_set_ptr());
	APolarHydrogenGraph fg(restype.graph(), boost::keep_all(), filter);
	return fg;
}

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

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
