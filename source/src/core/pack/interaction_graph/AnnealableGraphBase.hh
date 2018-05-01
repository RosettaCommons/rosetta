// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/InteractionGraphBase.hh
/// @brief  Base interface for annealable graphs.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_pack_interaction_graph_AnnealableGraphBase_hh
#define INCLUDED_core_pack_interaction_graph_AnnealableGraphBase_hh

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class AnnealableGraphBase : public utility::pointer::ReferenceCount
{
public:

	/// @brief Constructor.
	///
	AnnealableGraphBase();

	/// @brief Copy constructor.
	///
	AnnealableGraphBase( AnnealableGraphBase const &src );

	/// @brief Destructor.
	///
	virtual ~AnnealableGraphBase();

	// Graph property accessors.

	/// @brief Get the number of nodes in the graph.
	/// @details Must be implemented by derived classes.
	virtual int get_num_nodes() const = 0;

	/// @brief Get the number of states for a node in the graph.
	/// @details Must be implemented by derived classes.
	virtual int get_num_states_for_node(int node) const = 0;

	/// @brief Get the total number of states.
	/// @details Must be implemented by derived classes.
	virtual int get_num_total_states() const = 0;

	/// @brief Utility signal.
	/// @details Must be implemented by derived classes.
	virtual void prepare_for_simulated_annealing() = 0;

	/// @brief State initialization
	/// @details Must be implemented by derived classes.
	virtual void blanket_assign_state_0() = 0;

	/// @brief Is any state of any vertex unassigned?
	/// @details Must be implemented by derived classes.
	virtual bool any_vertex_state_unassigned() const = 0;

	/// @brief Explicit state modification for a node.
	/// @details Must be implemented by derived classes.
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state) = 0;

	/// @brief Explicit state modification for the network.
	/// @details Must be implemented by derived classes.
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states) = 0;

	/// @brief Consider switching node "node_ind" to state "new_state".
	/// @details Must be implemented by derived classes.
	/// @param[in] node_ind The node index.
	/// @param[in] new_state The state index that we are considering switching TO.
	/// @param[out] delta_energy The change in energy that results from the switch under consideration, computed by this function.
	/// @param[out] prev_energy_for_node The energy of this node prior to the substitutio, returned by this function.
	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	) = 0;

	/// @brief Accept a considered substitution (i.e. make the substitution permanent; commit it).
	/// @details Must be implemented by derived classes.
	virtual core::PackerEnergy commit_considered_substitution() = 0;

	/// @brief Get the energy resulting from the current set of state assignments.
	/// @details Must be implemented by derived classes.
	virtual core::PackerEnergy get_energy_current_state_assignment() = 0;

	/// @brief Set an error threshold.
	/// @details Must be implemented by derived classes.
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE ) = 0;

	/// @brief Provide the opportunity for an AnnealableGraph to clean up cached data in the pose or inside itself after packing.
	/// @details Base class function does nothing; may be overridden in derived classes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual void clean_up_after_packing( core::pose::Pose & pose );

};

}
}
}

#endif
