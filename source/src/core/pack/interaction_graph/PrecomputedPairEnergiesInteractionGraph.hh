// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh
/// @brief  Precomputed interaction graph class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_PrecomputedPairEnergiesInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_PrecomputedPairEnergiesInteractionGraph_hh

// Unit Headers
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>

// ObjexxFCL Headers

// Utility Headers

namespace core {
namespace pack {
namespace interaction_graph {

class PrecomputedPairEnergiesNode : public FixedBBNode
{
public:
	virtual ~PrecomputedPairEnergiesNode() {}

	PrecomputedPairEnergiesNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states)
	:
		FixedBBNode( owner, node_id, num_states )
	{}

};

class PrecomputedPairEnergiesEdge : public FixedBBEdge
{
public:
	virtual ~PrecomputedPairEnergiesEdge() {}

	PrecomputedPairEnergiesEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind)
	:
		FixedBBEdge( owner, first_node_ind, second_node_ind )
	{}

	virtual
	void add_to_two_body_energy(int const, int const, core::PackerEnergy const) = 0;

	virtual
	void add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array ) = 0;

	virtual
	void set_two_body_energy(int const, int const, core::PackerEnergy const) = 0;

	virtual
	void clear_two_body_energy(int const, int const) = 0;

};

class PrecomputedPairEnergiesInteractionGraph : public FixedBBInteractionGraph
{
public:
	virtual ~PrecomputedPairEnergiesInteractionGraph() {}
	PrecomputedPairEnergiesInteractionGraph( int num_nodes )
	:
		FixedBBInteractionGraph( num_nodes )
	{}


	/// @brief interface for PrecomputedPairEnergiesEdge::add_to_two_body_energies
	void add_to_two_body_energies_for_edge
	(
		int node1,
		int node2,
		ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array
	);

	/// @brief interface to PrecomputedPairEnergiesEdge::add_to_two_body_energies
	void add_to_two_body_energies_for_edge
	(
		int node1,
		int node2,
		int state_node1,
		int state_node2,
		core::PackerEnergy const two_body_energy
	);

	/// @brief interface to PDEdge::set_two_body_energy
	void set_two_body_energy_for_edge
	(
		int node1,
		int node2,
		int state_node1,
		int state_node2,
		core::PackerEnergy const two_body_energy
	);

	/// @brief interface to PDEdge::clear_two_body_energy
	void clear_two_body_energy_for_edge
	(
		int node1,
		int node2,
		int state_node1,
		int state_node2
	);

	virtual void declare_edge_energies_final(int node1, int node2);

};

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif
