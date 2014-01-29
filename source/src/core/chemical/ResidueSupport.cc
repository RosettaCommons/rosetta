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


//#ifndef INCLUDED_core_chemical_ResidueSupport_HH
//#define INCLUDED_core_chemical_ResidueSupport_HH

// Package Headers
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/graph/Graph.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>


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


/*	std::map<VD, std::map<VD, Size> > distance_matrix; //setup the matrix where Size is the distance between two vertex
	std::map<ED, Size> weight_map; //setup property (weight) map for the shortest distance


	for(EIterPair ep = boost::edges(graph_); ep.first != ep.second; ++ep.first){
		weight_map[*ep.first] = 1; //because the weight is 1. There is no distance asspciated with residues
	}

	//boost::weight_map(weight_map)
	boost::associative_property_map<std::map<ED, Size> > weight_property_map(weight_map);
	boost::floyd_warshall_all_pairs_shortest_paths(graph_, distance_matrix, boost::weight_map(weight_property_map) );

	Size const inf( 12345678 ); //assumption: fewer than 12 million nodes in the graph.
	//assert( num_nodes_ < inf );
	ObjexxFCL::FArray2D_int distance_table( natoms(), natoms(), inf);

	for(std::map<VD, std::map<VD, Size> >::const_iterator it = distance_matrix.begin(); it != distance_matrix.end(); ++it){
		VD vertex1 = it->first;
		std::map<VD, Size> value(it->second);
		for(std::map<VD, Size>::const_iterator second_it = value.begin(); second_it != value.end(); ++second_it){
			VD vertex2 = second_it->first;
			Size distance = second_it->second;
			distance_table(vd_to_index_.find(vertex1)->second , vd_to_index_.find(vertex2)->second ) = distance;
			distance_table(vd_to_index_.find(vertex2)->second, vd_to_index_.find(vertex1)->second) = distance;
		}
	}

    	std::cout << "Path Distances funciton call!" << std::endl;

	for(core::Size i = 1; i <= natoms(); ++i){
		for(core::Size ii =1; ii <= natoms(); ++ii){
			//std::cout << distance_table[i][ii];
			std::cout << name3() << " " << atom_name(i) << " " << atom_name(ii) << " " << distance_table(i, ii) << std::endl;
		}
	}

	return distance_table;*/


}

}
}

//#endif
