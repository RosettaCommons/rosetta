// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh
/// @brief  Annealable evaluator for score types evaluated over explicit list of residues.
/// @author Alex Ford (fordas@uw.edu)
//
#ifndef INCLUDED_core_pack_interaction_graph_ResidueArrayAnnealingEvaluator_hh
#define INCLUDED_core_pack_interaction_graph_ResidueArrayAnnealingEvaluator_hh

#include <utility>
#include <list>

// Project Headers
#include <core/types.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <utility/graph/Graph.fwd.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.fwd.hh>
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>

#include <map>

namespace core {
namespace pack {
namespace interaction_graph {

class ResidueArrayAnnealingEvaluator : public core::pack::interaction_graph::AnnealableGraphBase
{
public:
	/// @brief Constructor
	///
	ResidueArrayAnnealingEvaluator();

	/// @brief Destructor
	//.
	~ResidueArrayAnnealingEvaluator();

	/// @brief Copy constructor
	///
	ResidueArrayAnnealingEvaluator( ResidueArrayAnnealingEvaluator const &src );

	// @brief Initialization to be called by InterationGraphFactory during
	// AnneableGraph preparation.
	// @remark Initialize extracts all ResidueArrayAnnealingEvaluator-compatible score types
	// from the target score function and enables any with non-zero weights.
	void initialize(
		core::scoring::ScoreFunction const & score_function,
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets & rotamer_sets,
		utility::graph::GraphCOP
	);

	// Graph property accessors.

	/// @brief Get the number of nodes.
	/// @details This lies a little bit.  It only returns the number of independent nodes.
	virtual int get_num_nodes() const;

	/// @brief Get the number of states for a specific node.
	/// @param[in] node Index of the node.
	virtual int get_num_states_for_node(int node) const;

	/// @brief Get the total number of states for all nodes.
	/// @details This lies a little bit.  It only returns the total number of states for the independent nodes.
	virtual int get_num_total_states() const;

	/// @brief Utility signal.
	///
	virtual void prepare_for_simulated_annealing();

	/// @brief State initialization: set all nodes to state zero.
	///
	virtual void blanket_assign_state_0();

	/// @brief Are there any nodes unassigned?
	///
	virtual bool any_vertex_state_unassigned() const;

	/// @brief Explicit state modification: set a particular node to a particular state.
	/// @param[in] node_ind Index of the node to modify.
	/// @param[in] new_state Index of the state that we're setting this node TO.
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state);

	/// @brief Set states for all nodes across the network.
	/// @param[in] node_states Fortran-style 1-array of state indices for all nodes in the network.
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states);

	/// @brief Consider setting a particular node to a particular state.
	/// @param[in] node_ind Index of the node to be potentially modified.
	/// @param[in] new_state Index of the state that we might set this node to.
	/// @param[out] delta_energy The change in energy that would result from the substitution, returned by this function.
	/// @param[out] prev_energy_for_node The previous energy for this node, prior to the substitution, returned by this function.
	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	);

	/// @brief Set the node that we were considering to the state that we were considering (i.e. commit the change).
	///
	virtual core::PackerEnergy commit_considered_substitution();

	/// @brief Get the energy fo the current state.
	///
	virtual core::PackerEnergy get_energy_current_state_assignment();

	/// @brief Set error threshold.
	/// @param[in] deltaE Error threshold value to set.
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE );
	/// @brief Do any energy methods have a nonzero weight?
	///
	virtual bool has_methods() { return !weighted_energy_methods_.empty(); }

private:

	/// @brief Calculate the energy given the vector of residue owning pointers.
	/// @param[in] resvect 1-vector of const-owning pointers to Residue objects representing current state.
	core::Real calculate_weighted_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, int const substitution_position=0 );

	/// @brief Initialize symmetry information.
	/// @details Called by initialize().
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void initialize_symmetry_info( core::pose::Pose const &pose );

	/// @brief Sets the current consideration.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void set_consideration( int const node_ind, int const node_resid, int const new_state, utility::vector1< std::pair< int, core::conformation::ResidueCOP> > & unset_info );

	/// @brief Clears the current consideration.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void unset_consideration( utility::vector1< std::pair < int, core::conformation::ResidueCOP> > const & unset_info );


	/// @brief Abandon the substitution that was being considered.
	///
	void clear_consideration();

private:

	typedef std::pair< core::PackerEnergy, core::scoring::annealing::ResidueArrayAnnealableEnergyOP > WeightedMethodPair;

	/// @brief List of energy methods and their weights.
	///
	std::list< WeightedMethodPair > weighted_energy_methods_;

	/// @brief Const-owning pointers to the residues of the source pose.
	///
	utility::vector1< core::conformation::ResidueCOP > source_pose_residues_;

	/// @brief Const-owning pointers to the current residues.
	///
	utility::vector1< core::conformation::ResidueCOP > current_residues_;

	/// @brief Current energy.
	///
	core::PackerEnergy current_energy_;

	/// @brief The node being considered for alteration.
	///
	int considered_node_;

	/// @brief The state that we're considering changing this node TO.
	///
	int considered_state_;

	/// @brief The energy that would result from the considered substitution.
	///
	core::PackerEnergy considered_energy_;

	/// @brief The rotamer sets for each node (vector of const-owning pointers to RotamerSet objects.
	///
	utility::vector1< core::pack::rotamer_set::RotamerSetCOP > per_node_rotamer_sets_;

	/// @brief Is this a symmetric packing job?
	///
	bool is_symmetric_;

	/// @brief Is this a symmetric packing job with mirror symmetry?
	///
	bool has_mirror_symm_;

	/// @brief Number of independent nodes (not symmetry copies).
	///
	core::Size num_indep_nodes_;

	/// @brief Map of dependent nodes.
	/// @details Maps controlling node index -> vector of dependent node indices.
	std::map < int, utility::vector1< core::Size > > dependent_node_map_;

	/// @brief Map of dependent residues.
	/// @details Maps controlling node index -> vector of dependent residue indices.
	std::map < int, utility::vector1< core::Size > > dependent_residue_map_;

};

}
}
}

#endif
