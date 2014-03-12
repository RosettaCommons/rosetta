// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/ResidueSupport.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley

// Package Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
// Project Headers
#include <core/graph/Graph.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <utility/graph/RingDetection.hh>

// ObjexxFCL Headers
// Commented by inclean daemon #include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace chemical {

ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res )
{
	//return res.get_residue_path_distances();
	using namespace graph;
	Graph g;

	g.set_num_nodes( res.natoms() );
	for ( uint ii = 1; ii <= res.natoms(); ++ii )
	{
		AtomIndices const ii_bonded = res.nbrs( ii  );
		for ( Size jj = 1; jj <= ii_bonded.size(); ++jj)
		{
			if ( ii_bonded[ jj ] > ii )
				g.add_edge( ii, ii_bonded[ jj ] );
		}
	}
	return g.all_pairs_shortest_paths();

}


LightWeightResidueGraph convert_residuetype_to_light_graph(ResidueType const & res){
	//this is a const reference because vertices change when  you copy the graph
	const  core::chemical::ResidueGraph & full_residue_graph = res.graph(); //get the boost graph structure from residuetype

	LightWeightResidueGraph lwrg;

	//set up the mapping between VD and ED from ResidueGraph for LightWeightResidueGraph
	boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);

	//map between the LightWeightResidueGraph vertex and the ResidueGraph vertex. Used when adding edges
	std::map<core::chemical::VD, lwrg_VD> map_of_vertex;
	//first, add all the vertex to light weight residue graph, and map the property maps to point at the ResidueGraph
	for(core::chemical::VIterPair vp = boost::vertices(full_residue_graph); vp.first != vp.second; ++vp.first){
		VD const & full_residue_graph_vd = *vp.first;
		lwrg_VD lwrg_vd = boost::add_vertex(lwrg);
		lwrg_vd_to_VD[lwrg_vd] = full_residue_graph_vd; //set property maps
		map_of_vertex[full_residue_graph_vd] = lwrg_vd; //set mapping for the edges
	}

	//now we add the edges between the vertex
	for(EIterPair ep = boost::edges(full_residue_graph); ep.first != ep.second; ++ep.first){
		VD source = boost::source(*ep.first, full_residue_graph);
		VD target = boost::target(*ep.first, full_residue_graph);
		lwrg_VD source_lwrg_VD = map_of_vertex[source];
		lwrg_VD target_lwrg_VD = map_of_vertex[target];
		lwrg_ED e_added;
		bool added;
		boost::tie(e_added, added) = boost::add_edge(source_lwrg_VD, target_lwrg_VD,lwrg );
		if(added){ //only add bonds once!
			lwrg_ed_to_ED[e_added] = *ep.first;
		}
	}
	assert(boost::num_vertices(lwrg) == full_residue_graph.num_vertices()); //fail if the number of vertex are not the same
	assert(boost::num_edges(lwrg) == full_residue_graph.num_edges()); //fail if the number of edges are not the same


	//boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	//boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);


	utility::graph::RingDetection ring_detect(lwrg); //initialize the ring detector. Automatically assigns rings
	utility::vector1<utility::vector1<lwrg_VD> > rings = ring_detect.GetRings(); //these are the path of the rings


	//our light weight graph has been made. Now return it!
	return lwrg;
}


}
}

//#endif
