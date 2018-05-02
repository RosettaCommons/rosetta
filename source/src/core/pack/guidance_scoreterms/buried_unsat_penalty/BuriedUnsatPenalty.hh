// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.hh
/// @brief Headers for an EnergyMethod that gives a penalty for buried unsatisfied hydrogen bond donors and acceptors.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_BuriedUnsatPenalty_hh
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_BuriedUnsatPenalty_hh

// Unit headers
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.fwd.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.fwd.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.fwd.hh>

// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

// Utility headers

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>

// Forward declarations for friendship:
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltyTests.fwd.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {

// Uncomment the following line to do a full recomputation at every step, to confirm correct decrementation.  This is SLOW and should
// not be used for production runs:
//#define SLOW_BUNSAT_PENALTY_ACCURACY_CHECKS_ENABLED

/// @brief BuriedUnsatPenalty, an EnergyMethod that gives a penalty for buried unsatisfied hydrogen bond donors or acceptors.
/// @details This class is derived from base class WholeStructureEnergy, which is meaningful only on entire structures.
/// These EnergyMethods do all of their work in the "finalize_total_energy" section of scorefunction evaluation.
class BuriedUnsatPenalty : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {

	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	BuriedUnsatPenalty( core::scoring::methods::EnergyMethodOptions const &options );

	/// @brief Default destructor.
	///
	virtual ~BuriedUnsatPenalty();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief BuriedUnsatPenalty is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief BuriedUnsatPenalty is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.  The update_residue_neighbors() function of the pose
	/// must be called first.
	void finalize_total_energy( core::pose::Pose & pose, core::scoring::ScoreFunction const & sfxn, core::scoring::EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	/// @note Outside of the context of the packer, this doesn't behave as expected (and finalize_total_energy() should
	/// be called instead).
	core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, core::Size const substitution_position = 0 ) const override;

	/// @brief What to do when a substitution that was considered is accepted.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void commit_considered_substitution() override;

	/// @brief Get a summary of all loaded data.
	///
	void report() const;

	/// @brief Cache data from the pose in this EnergyMethod in anticipation of packing.
	///
	void set_up_residuearrayannealableenergy_for_packing( core::pose::Pose &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn) override;

	/// @brief Delete cached data from the pose and from this EnergyMethod after packing.
	///
	void clean_up_residuearrayannealableenergy_after_packing( core::pose::Pose &pose ) override;

	/// @brief Disable this scoreterm during minimization trajectory.
	void setup_for_minimizing( core::pose::Pose & pose, core::scoring::ScoreFunction const & sfxn, core::kinematics::MinimizerMapBase const & minmap) const override;

	/// @brief Re-enable this scoreterm after a minimization trajectory.
	void finalize_after_minimizing( core::pose::Pose & pose) const override;

public:

	/******************************************************
	Public member functions specific to this energy method.
	******************************************************/

	/// @brief Given the counts of various unsaturateds, return a penalty value.
	core::Real compute_penalty( core::Size const unsat_acceptor_count, core::Size const unsat_donor_count, core::Size const unsat_acceptor_and_donor_count, core::Size const oversat_acceptor_count, core::Size const oversat_donor_count, core::Size const oversat_acceptor_and_donor_count ) const;

	/// @brief Given a residue vector, calculate the penalty energy using a reference graph.
	core::Real calculate_penalty_once_almost_from_scratch_using_reference_graph( utility::vector1< core::conformation::ResidueCOP > const & resvect, graph::BuriedUnsatPenaltyGraph const &reference_graph, core::Size const symm_multiplier ) const;

	/// @brief Given a pose, calculate the penalty energy.
	core::Real calculate_penalty_once_from_scratch( core::pose::Pose const &pose ) const;

	/// @brief Provide Pymol commands to colour the pose grey, non-buried donor and acceptor groups cyan, and buried acceptor
	/// and donor groups orange.  Useful for debugging degree of burial.
	/// @details To use, pass in a pose.  If this graph contains residues corresponding to those in the pose, commands for colouring
	/// them will be written out.  Calls BuriedUnsatPenaltyGraph::provide_pymol_commands_to_show_groups().  Must be called only after
	/// set_up_residuearrayannealableenergy_for_packing.
	void provide_pymol_commands_to_show_groups( std::ostream &out, core::pose::Pose const &pose ) const;

private:

	/******************
	Private functions:
	******************/

	/// @brief Called from set_up_residuearrayannelableenergy_for_packing().  Initializes the graph structure representing the current state
	/// during the packing trajectory.
	void initialize_curstate_graph( core::pose::Pose const & pose );

	/// @brief Initialize the graph structure representing the current state during the packing trajectory from a vector of residues.  Called from
	/// the first step of calculate_energy().
	void initialize_curstate_graph( utility::vector1< core::conformation::ResidueCOP > const &resvect );

	/// @brief Given a list of changed node indices and a graph of the current state, determine which nodes share edges with the changed nodes, and add
	/// their indices to a list of partners.
	/// @details Don't add indices that are in the changed_node_indices list or already in the changed_node_partners list.
	void add_to_list_of_partners_of_changed_nodes( graph::BuriedUnsatPenaltyGraphCOP curstate_graph, utility::vector1< core::Size > const & changed_node_indices, utility::vector1< core::Size > & changed_node_partners ) const;

	/// @brief Increment the counts based on the current state of the curstate_graph_ and the current changed_node_indices_ and changed_node_partners_ vectors.
	void increment_counts(
		core::Size & unsat_acceptor_count_lastconsidered,
		core::Size & unsat_donor_count_lastconsidered,
		core::Size & unsat_acceptor_and_donor_count_lastconsidered,
		core::Size & oversat_acceptor_count_lastconsidered,
		core::Size & oversat_donor_count_lastconsidered,
		core::Size & oversat_acceptor_and_donor_count_lastconsidered,
		utility::vector1< core::Size > const & changed_node_indices,
		utility::vector1< core::Size > const & changed_node_partners
	) const;

	/// @brief Decrement the counts based on the current state of the curstate_graph_ and the current changed_node_indices_ and changed_node_partners_ vectors.
	void decrement_counts(
		core::Size & unsat_acceptor_count_lastconsidered,
		core::Size & unsat_donor_count_lastconsidered,
		core::Size & unsat_acceptor_and_donor_count_lastconsidered,
		core::Size & oversat_acceptor_count_lastconsidered,
		core::Size & oversat_donor_count_lastconsidered,
		core::Size & oversat_acceptor_and_donor_count_lastconsidered,
		utility::vector1< graph::BuriedUnsatPenaltyNodeDataCOP > const &old_data,
		utility::vector1< core::Size > const & changed_node_indices,
		utility::vector1< core::Size > const & changed_node_partners
	) const;

	/******************
	Private variables:
	******************/

	/// @brief Is this scoreterm currently disabled (e.g. for a minimization trajectory)?
	mutable bool disabled_;

	/// @brief Pointer to buried unsaturated penalty graph in the pose.  Stored during packing
	/// to avoid repeated lookups.
	graph::BuriedUnsatPenaltyGraphCOP unsat_graph_;

	/// @brief Options for the BruiedUnsatPenaltyGraph.
	graph::BuriedUnsatPenaltyGraphOptionsOP graph_options_;

	/// @brief Options for hydrogen bonds.
	core::scoring::hbonds::HBondOptionsOP hbond_options_;

	/// @brief Is this pose symmetric?  Used during packing only.
	bool symmetric_;

	/// @brief How many residues are in the asymmetric unit of a symmetric pose (symmetric case) or in a pose (asymmetric case)?
	/// Used during packing only.
	core::Size nres_;

	/// @brief The number of symmetry copies.  Used for packing only.
	core::Size num_symmetric_copies_;

	/******************
	Private variables used in packing:
	******************/

	/// @brief The graph representing the last accepted state during the packer trajectory.
	mutable graph::BuriedUnsatPenaltyGraphOP curstate_graph_;

	/// @brief The unsat acceptor counts last accepted.
	mutable core::Size unsat_acceptor_count_lastaccepted_;

	/// @brief The unsat donor counts last accepted.
	mutable core::Size unsat_donor_count_lastaccepted_;

	/// @brief The unsat acceptor and donor counts last accepted.
	mutable core::Size unsat_acceptor_and_donor_count_lastaccepted_;

	/// @brief The oversat acceptor counts last accepted.
	mutable core::Size oversat_acceptor_count_lastaccepted_;

	/// @brief The oversat donor counts last accepted.
	mutable core::Size oversat_donor_count_lastaccepted_;

	/// @brief The oversat acceptor and donor counts last accepted.
	mutable core::Size oversat_acceptor_and_donor_count_lastaccepted_;

#ifdef SLOW_BUNSAT_PENALTY_ACCURACY_CHECKS_ENABLED
	/// @brief A copy of the pose, used for slow accuracy checks only.
	core::pose::PoseOP temppose_;

	/// @brief A counter for the step in the packer trajectory
	mutable core::Size packer_step_;
#endif

}; //BuriedUnsatPenalty class

} // buried_unsat_penalty
} // guidance_scoreterms
} // pack
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
