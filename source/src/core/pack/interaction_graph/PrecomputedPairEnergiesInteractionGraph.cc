// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.cc
/// @brief  Precomputed interaction graph class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add support for multithreaded interaction graph precomputation.

// Associated Headers
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>

// Core Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// STL Headers
#include <iostream>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#include <basic/thread_manager/RosettaThreadManager.hh>

using namespace ObjexxFCL;

static basic::Tracer TR( "core.pack.interaction_graph.PrecomputedPairEnergiesInteractionGraph", basic::t_info );

namespace core {
namespace pack {
namespace interaction_graph {


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param oversized_res_res_energy_array - [in] - the large table
///  of rotamer pair energies
/// @param [in] use_threadsafe_method If true, we do NOT cache the edge pointer for faster
/// future lookups.
void PrecomputedPairEnergiesInteractionGraph::add_to_two_body_energies_for_edge
(
	int node1,
	int node2,
	FArray2< core::PackerEnergy > const & oversized_res_res_energy_array,
	bool const use_threadsafe_method /*= false*/
)
{
	auto* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2, use_threadsafe_method );
	if ( edge == nullptr ) {
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
	auto* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if ( edge == nullptr ) {
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
	auto* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if ( edge == nullptr ) {
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
	auto* edge =
		(PrecomputedPairEnergiesEdge*) find_edge( node1, node2 );
	if ( edge == nullptr ) return;
	edge->set_two_body_energy( state_node1, state_node2, 0 );
}


/// @brief Given two nodes, compute the energy of all interacting rotamers and add those
/// energies to the edge between the nodes.  This is suitable for use in a multithreaded
/// context.
/// @details The edge must already exist.  The "do_orient" and "symmetric_swap" inputs are
/// only used in the symmetric case, and only control the generation of temporary, oriented
/// rotamersets from rotsets.
/// @note If this is a symmetric packing case, symminfo must be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PrecomputedPairEnergiesInteractionGraph::set_twobody_energies_multithreaded(
	std::pair< core::Size, core::Size > const &nodes,
	std::pair< core::pack::rotamer_set::RotamerSetCOP, core::pack::rotamer_set::RotamerSetCOP > rotsets,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
#ifdef MULTI_THREADED
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_info,
#else
	basic::thread_manager::RosettaThreadAssignmentInfo const &,
#endif
	bool const finalize_edges,
	core::conformation::symmetry::SymmetryInfoCOP symminfo,
	std::pair< core::Size, core::Size > const &resids,
	bool const do_orient,
	bool const symmetric_swap
) {
#ifdef MULTI_THREADED
	if ( TR.Debug.visible() ) {
		TR.Debug << "Thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " (of " << thread_info.get_assigned_total_thread_count() << " threads assigned to interaction graph precomputation) computing" << (symminfo == nullptr ? " " : " symmetric ") << "interaction between nodes " << nodes.first  << " and " << nodes.second << "." << std::endl;
	}
#endif //MULTI_THREADED

	// If we have a intersubunit interaction, then calculate the interaction
	// here instead of the intrasubunit interaction. If jj == 0 then we will calculate
	// intrasubunit interaction. If we swapped the order we have to orient ii instead of jj
	if ( symminfo != nullptr && do_orient ) {
		if ( symmetric_swap ) {
			core::pack::rotamer_set::RotamerSetOP rotated_ii_rotset(
				core::pack::rotamer_set::symmetry::SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(pose,rotsets.first,resids.second) );
			sfxn.prepare_rotamers_for_packing( pose, *rotated_ii_rotset );
			rotsets.first = rotated_ii_rotset;
		} else {
			core::pack::rotamer_set::RotamerSetOP rotated_jj_rotset(
				core::pack::rotamer_set::symmetry::SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(pose,rotsets.second,resids.second) );
			sfxn.prepare_rotamers_for_packing( pose, *rotated_jj_rotset );
			rotsets.second = rotated_jj_rotset;
		}
	}

	FArray2D< core::PackerEnergy > pair_energy_table( rotsets.second->num_rotamers(), rotsets.first->num_rotamers(), 0.0 );
	sfxn.evaluate_rotamer_pair_energies( *(rotsets.first), *(rotsets.second), pose, pair_energy_table );

	if ( symminfo != nullptr ) {
		pair_energy_table *= symminfo->score_multiply(resids.first, resids.second);
	}

	add_to_two_body_energies_for_edge( nodes.first, nodes.second, pair_energy_table, true );

	if ( finalize_edges && ! sfxn.any_lr_residue_pair_energy(pose, nodes.first, nodes.second) ) {
		declare_edge_energies_final( nodes.first, nodes.second, true );
	}
}

/// @brief Given two nodes and a long-range energy, compute the long-range energy of all
/// interacting rotamers and add those energies to the edge between the nodes.  This is
/// suitable for use in a multithreaded context.
/// @details The edge must already exist.  The "do_orient" and "symmetric_swap" inputs are
/// only used in the symmetric case, and only control the generation of temporary, oriented
/// rotamersets from rotsets.
/// @note If this is a symmetric packing case, symminfo must be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PrecomputedPairEnergiesInteractionGraph::add_longrange_twobody_energies_multithreaded(
	std::pair< core::Size, core::Size > const &nodes,
	std::pair< core::pack::rotamer_set::RotamerSetCOP, core::pack::rotamer_set::RotamerSetCOP > rotsets,
	core::pose::Pose const & pose,
	core::scoring::methods::LongRangeTwoBodyEnergyCOP lr_energy,
	core::scoring::ScoreFunction const & sfxn,
#ifdef MULTI_THREADED
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_info,
#else
	basic::thread_manager::RosettaThreadAssignmentInfo const &,
#endif
	bool const finalize_edges,
	core::conformation::symmetry::SymmetryInfoCOP symminfo,
	std::pair< core::Size, core::Size > const &resids,
	bool const do_orient,
	bool const symmetric_swap
) {
#ifdef MULTI_THREADED
	if ( TR.Debug.visible() ) {
		TR.Debug << "Thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " (of " << thread_info.get_assigned_total_thread_count() << " threads assigned to interaction graph precomputation) computing a" << (symminfo == nullptr ? " " : " symmetric ") << "long-range energy between nodes " << nodes.first  << " and " << nodes.second << "." << std::endl;
	}
#endif //MULTI_THREADED

	// we have a intersubunit interaction then calculate the interaction
	// here instead of the intrasubunit interaction. If jj == 0 then we will calculate
	// intrasubunit interaction. If we swapped the order we have to orient ii instead of jj
	if ( symminfo != nullptr && do_orient ) {
		if ( symmetric_swap ) {
			core::pack::rotamer_set::RotamerSetOP rotated_ii_rotset(
				core::pack::rotamer_set::symmetry::SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(pose, rotsets.first , resids.second) );
			sfxn.prepare_rotamers_for_packing( pose, *rotated_ii_rotset );
			rotsets.first = rotated_ii_rotset;
		} else {
			core::pack::rotamer_set::RotamerSetOP rotated_jj_rotset(
				core::pack::rotamer_set::symmetry::SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(pose, rotsets.second , resids.second) );
			sfxn.prepare_rotamers_for_packing( pose, *rotated_jj_rotset );
			rotsets.second = rotated_jj_rotset;
		}
	}

	FArray2D< core::PackerEnergy > pair_energy_table( rotsets.second->num_rotamers(), rotsets.first->num_rotamers(), 0.0 );
	lr_energy->evaluate_rotamer_pair_energies( *(rotsets.first), *(rotsets.second), pose, sfxn, sfxn.weights(), pair_energy_table );

	if ( symminfo != nullptr ) {
		pair_energy_table *= symminfo->score_multiply(resids.first, resids.second);
	}

	add_to_two_body_energies_for_edge( nodes.first, nodes.second, pair_energy_table, true );
	if ( finalize_edges ) declare_edge_energies_final( nodes.first, nodes.second, true );
}

/// @brief call this if you're done storing energies in an edge - it will reduce
/// the memory usage for that edge if possible
///
/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param use_threadsafe_lookup - [in] - If true, we don't cache the last edge
/// found for faster future lookups.  False by default.
void PrecomputedPairEnergiesInteractionGraph::declare_edge_energies_final
(
	int node1,
	int node2,
	bool const use_threadsafe_lookup
)
{
	auto* edge = static_cast<PrecomputedPairEnergiesEdge*>( find_edge( node1, node2, use_threadsafe_lookup ) );
	if ( edge == nullptr ) {
		return;
	}
	edge->declare_energies_final();
}

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core
