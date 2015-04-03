// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.cc
/// @brief  Precomputed interaction graph class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Rosetta Headers
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// STL Headers
#include <iostream>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param oversized_res_res_energy_array - [in] - the large table
/// 	of rotamer pair energies
void PrecomputedPairEnergiesInteractionGraph::add_to_two_body_energies_for_edge
(
	int node1,
	int node2,
	FArray2< core::PackerEnergy > const & oversized_res_res_energy_array
)
{
	PrecomputedPairEnergiesEdge* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if (edge == NULL){
		std::cerr << "WARNING:: you've input edge energies for an edge that does not exist" << std::endl;
		return;
	}
	edge->add_to_two_body_energies( oversized_res_res_energy_array );
}


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param state_node1 - [in] - state on smaller-indexed node
/// @param state_node2 - [in] - state on larger-indexed node
/// @param two_body_energy  - [in] - the energy for this state pair
void
PrecomputedPairEnergiesInteractionGraph::add_to_two_body_energies_for_edge
(
	int node1,
	int node2,
	int state_node1,
	int state_node2,
	core::PackerEnergy const two_body_energy
)
{
	PrecomputedPairEnergiesEdge* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if (edge == NULL) {
		std::cerr << "WARNING:: you've input edge energies for an edge that does not exist" << std::endl;
		return;
	}
	edge->add_to_two_body_energy( state_node1, state_node2, two_body_energy );
}


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param state_node1 - [in] - state on smaller-indexed node
/// @param state_node2 - [in] - state on larger-indexed node
/// @param two_body_energy  - [in] - the energy for this state pair
void PrecomputedPairEnergiesInteractionGraph::set_two_body_energy_for_edge(
	int node1,
	int node2,
	int state_node1,
	int state_node2,
	core::PackerEnergy const two_body_energy
)
{
	PrecomputedPairEnergiesEdge* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if (edge == NULL) {
		std::cerr << "WARNING:: you've input edge energies for an edge that does not exist" << std::endl;
		return;
	}
	edge->set_two_body_energy( state_node1, state_node2, two_body_energy );
}

/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param state_node1 - [in] - state on smaller-indexed node
/// @param state_node2 - [in] - state on larger-indexed node
///
/// @author jk
///
void PrecomputedPairEnergiesInteractionGraph::clear_two_body_energy_for_edge(
	int node1,
	int node2,
	int state_node1,
	int state_node2
)
{
	PrecomputedPairEnergiesEdge* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if (edge == NULL) return;
	edge->set_two_body_energy( state_node1, state_node2, 0 );
}


/// @brief call this if you're done storing energies in an edge - it will reduce
/// the memory usage for that edge if possible
///
/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
void PrecomputedPairEnergiesInteractionGraph::declare_edge_energies_final
(
	int node1,
	int node2
)
{
	PrecomputedPairEnergiesEdge* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if (edge == NULL) {
		return;
	}
	edge->declare_energies_final();
}

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core
