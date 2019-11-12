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
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add support for multithreaded interaction graph precomputation.


#ifndef INCLUDED_core_pack_interaction_graph_PrecomputedPairEnergiesInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_PrecomputedPairEnergiesInteractionGraph_hh

// Unit Headers
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>

// Core Headers
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>

// ObjexxFCL Headers

// Utility Headers

// Basic Headers
#include <basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class PrecomputedPairEnergiesNode : public FixedBBNode
{
public:
	~PrecomputedPairEnergiesNode() override {}

	PrecomputedPairEnergiesNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states)
	:
		FixedBBNode( owner, node_id, num_states )
	{}

public:

	/// @brief Get the one-body energy for a state of this node.
	/// @details Must be implmented by derived classes.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	virtual core::PackerEnergy get_one_body_energy( int state ) const = 0;

};

class PrecomputedPairEnergiesEdge : public FixedBBEdge
{
public:
	~PrecomputedPairEnergiesEdge() override {}

	PrecomputedPairEnergiesEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind)
	:
		FixedBBEdge( owner, first_node_ind, second_node_ind )
	{}


	/// @brief Get the two-body energy for two states of the nodes connected to this edge.
	/// @details Must be implmented by derived classes.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	core::PackerEnergy get_two_body_energy( int const first_state, int const second_state ) const override = 0;

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
	~PrecomputedPairEnergiesInteractionGraph() override {}
	PrecomputedPairEnergiesInteractionGraph( int num_nodes )
	:
		FixedBBInteractionGraph( num_nodes )
	{}


	/// @brief interface for PrecomputedPairEnergiesEdge::add_to_two_body_energies
	/// @note If use_threadsafe_method is true, we don't cache the edge index for faster future
	/// lookups.
	void add_to_two_body_energies_for_edge
	(
		int node1,
		int node2,
		ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array,
		bool const use_threadsafe_method = false
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

	/// @brief Given two nodes, compute the energy of all interacting rotamers and add those
	/// energies to the edge between the nodes.  This is suitable for use in a multithreaded
	/// context.
	/// @details The edge must already exist.  The "do_orient" and "symmetric_swap" inputs are
	/// only used in the symmetric case, and only control the generation of temporary, oriented
	/// rotamersets from rotsets.
	/// @note If this is a symmetric packing case, symminfo must be provided.  The resid1 and
	/// resid2 values need only be provided in the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	set_twobody_energies_multithreaded(
		std::pair< core::Size, core::Size > const &nodes,
		std::pair< core::pack::rotamer_set::RotamerSetCOP, core::pack::rotamer_set::RotamerSetCOP > rotsets,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		basic::thread_manager::RosettaThreadAssignmentInfo const & thread_info,
		bool const finalize_edges,
		core::conformation::symmetry::SymmetryInfoCOP symminfo,
		std::pair< core::Size, core::Size > const &resids,
		bool const do_orient,
		bool const symmetric_swap
	);

	/// @brief Given two nodes and a long-range energy, compute the long-range energy of all
	/// interacting rotamers and add those energies to the edge between the nodes.  This is
	/// suitable for use in a multithreaded context.
	/// @details The edge must already exist.  The "do_orient" and "symmetric_swap" inputs are
	/// only used in the symmetric case, and only control the generation of temporary, oriented
	/// rotamersets from rotsets.
	/// @note If this is a symmetric packing case, symminfo must be provided.  The resid1 and
	/// resid2 values need only be provided in the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	add_longrange_twobody_energies_multithreaded(
		std::pair< core::Size, core::Size > const &nodes,
		std::pair< core::pack::rotamer_set::RotamerSetCOP, core::pack::rotamer_set::RotamerSetCOP > rotsets,
		core::pose::Pose const & pose,
		core::scoring::methods::LongRangeTwoBodyEnergyCOP lr_energy,
		core::scoring::ScoreFunction const & sfxn,
		basic::thread_manager::RosettaThreadAssignmentInfo const & thread_info,
		bool const finalize_edges,
		core::conformation::symmetry::SymmetryInfoCOP symminfo,
		std::pair< core::Size, core::Size > const &resids,
		bool const do_orient,
		bool const symmetric_swap
	);


	/// @brief call this if you're done storing energies in an edge - it will reduce
	/// the memory usage for that edge if possible
	///
	/// @param node1 - [in] - the index of the smaller-indexed node
	/// @param node2 - [in] - the index of the larger-indexed node
	/// @param use_threadsafe_lookup - [in] - If true, we don't cache the last edge
	/// found for faster future lookups.  False by default.
	virtual
	void
	declare_edge_energies_final(
		int node1,
		int node2,
		bool const use_threadsafe_lookup = false
	);


};

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif
