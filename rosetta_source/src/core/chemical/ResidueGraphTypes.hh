// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_ResiduePredicates_hh
#define INCLUDED_core_chemical_ResiduePredicates_hh

// Unit headers
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <core/chemical/AtomTypeSet.fwd.hh>

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

typedef boost::graph_traits<ResidueGraph>::vertex_iterator VIter;
typedef boost::graph_traits<ResidueGraph>::edge_iterator EIter;
typedef std::pair<VIter, VIter> VIterPair;

typedef boost::graph_traits<ResidueGraph>::out_edge_iterator OutEdgeIter;
//typedef boost::graph_traits<ResidueGraph>::in_edge_iterator InEdgeIter; // Out and in edges are the same in an undirected_graph
typedef std::pair<OutEdgeIter, OutEdgeIter> OutEdgeIterPair;

typedef std::map< std::string, VD > NameVDMap;
typedef std::pair<std::string, VD> NameVDPair;
typedef std::pair<NameVDMap::iterator, bool> NameVDInserted;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////////// PREDICATES for FILTERED GRAPHS ///////////////////
////////////////////////////////////////////////////////////

class HeavyAtomFilter{
public:
	HeavyAtomFilter(){};
	HeavyAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(atom_types){};
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

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

class HeavyAtomWithPolarHydrogensFilter{
public:
	HeavyAtomWithPolarHydrogensFilter(){};
	HeavyAtomWithPolarHydrogensFilter(HeavyAtomGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(atom_types){};
	bool operator()(VD const vd) const;
private:
	HeavyAtomGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<HeavyAtomGraph, boost::keep_all, HeavyAtomWithPolarHydrogensFilter> HeavyAtomWithPolarHydrogensGraph;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::vertex_iterator HeavyAtomWithPolarHydrogensVIter;
typedef boost::graph_traits<HeavyAtomWithPolarHydrogensGraph>::edge_iterator HeavyAtomWithPolarHydrogensEIter;

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

class HeavyAtomAcceptorFilter{
public:
	HeavyAtomAcceptorFilter(){};
	HeavyAtomAcceptorFilter(HeavyAtomGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(atom_types){};
	bool operator()(VD const vd) const;
private:
	HeavyAtomGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<HeavyAtomGraph, boost::keep_all, HeavyAtomAcceptorFilter> HeavyAtomAcceptorGraph;
typedef boost::graph_traits<HeavyAtomAcceptorGraph>::vertex_iterator HeavyAtomAcceptorVIter;
typedef boost::graph_traits<HeavyAtomAcceptorGraph>::edge_iterator HeavyAtomAcceptorEIter;

////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//class PolarHydrogenFilter{
//public:
//	PolarHydrogenFilter(){};
//	PolarHydrogenFilter(HeavyAtomGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(atom_types){};
//	bool operator()(VD const vd) const;
//private:
//	HeavyAtomGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
//	AtomTypeSetCAP atom_types_;
//};
//typedef boost::filtered_graph<HeavyAtomGraph, boost::keep_all, HeavyAtomAcceptorFilter> HeavyAtomAcceptorGraph;
//typedef boost::graph_traits<HeavyAtomAcceptorGraph>::vertex_iterator HeavyAtomAcceptorVIter;
//typedef boost::graph_traits<HeavyAtomAcceptorGraph>::edge_iterator HeavyAtomAcceptorEIter;

}
}
///////////////////////////////////////////////////////////////


#endif
