// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/MultiplexedAnnealableGraph.hh
/// @brief  Multiplexing packing graph container.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_pack_interaction_graph_MultiplexedAnnealableGraph_hh
#define INCLUDED_core_pack_interaction_graph_MultiplexedAnnealableGraph_hh

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/pack/interaction_graph/MultiplexedAnnealableGraph.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <list>

namespace core {
namespace pack {
namespace interaction_graph {

class MultiplexedAnnealableGraph : public AnnealableGraphBase
{
	typedef std::list< AnnealableGraphBaseOP > SubgraphContainer;
	typedef SubgraphContainer::iterator       iterator;
	typedef SubgraphContainer::const_iterator const_iterator;

public:
	/// @brief Default constructor.
	///
	MultiplexedAnnealableGraph();

	/// @brief Copy constructor.
	///
	MultiplexedAnnealableGraph(MultiplexedAnnealableGraph const & other);

	/// @brief Setup constructor.
	///
	MultiplexedAnnealableGraph(SubgraphContainer const & target_subgraphs);

	/// @brief Destructor.
	///
	virtual ~MultiplexedAnnealableGraph();

	// AnnealableGraphBase Implementation
public:
	// Graph property accessors.

	/// @brief Get the number of nodes in the graph.
	///
	int get_num_nodes() const override;

	/// @brief Get the number of states for a given node.
	/// @param[in] node The index of the node.
	int get_num_states_for_node(int node) const override;

	/// @brief Get the total number of states for all nodes in the graph.
	///
	int get_num_total_states() const override;

	/// @brief Utility signal.
	///
	void prepare_for_simulated_annealing() override;

	/// @brief State initialization.  Set all nodes to state zero.
	///
	void blanket_assign_state_0() override;

	/// @brief Are any states unassigned?
	///
	bool any_vertex_state_unassigned() const override;

	/// @brief Explicit state modification: set node "node_ind" to state "new_state".
	/// @param[in] node_ind The index of the node.
	/// @param[in] new_state The index of the state to which we're setting the node.
	core::PackerEnergy set_state_for_node(int node_ind, int new_state) override;

	/// @brief Set the states for the entire network.
	/// @param[in] node_states A Fortran-style 1-array of state indices for all nodes in the network.
	core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states) override;


	/// @brief Consider a change in state at a particular node.
	/// @param[in] node_ind The node index.
	/// @param[in] new_state The index of the state to which we're considering setting the node.
	/// @param[out] delta_energy The computed change in energy that would result from the substitution.
	/// @param[out] prev_energy_for_state The energy prior to the substituion, returned by this function.
	void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	) override;

	/// @brief Accept the considered substitution.
	///
	core::PackerEnergy commit_considered_substitution() override;

	/// @brief Get the total energy from the current states of the nodes in the network.
	///
	core::PackerEnergy get_energy_current_state_assignment() override;

	/// @brief Set error threshold.
	/// @param[in] deltaE Error threshold value to set.
	void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE ) override;

	/// @brief Provide the opportunity for an AnnealableGraph to clean up cached data in the pose or inside itself after packing.
	/// @details This version calls contained AnnealableGraph clean_up_after_packing() methods.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void clean_up_after_packing( core::pose::Pose & pose ) override;

	/// @brief Container accessor
	///
	SubgraphContainer subgraphs;
};

}
}
}

#endif
