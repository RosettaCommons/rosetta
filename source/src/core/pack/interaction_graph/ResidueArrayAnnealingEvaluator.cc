// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh
/// @brief  Annealable interface for score terms evaluated over an array of residues.
/// @author Alex Ford (fordas@uw.edu)
//

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/assert.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/graph/Graph.fwd.hh>

#include <boost/foreach.hpp>
#define foreach_         BOOST_FOREACH
#include <utility>

namespace core {
namespace pack {
namespace interaction_graph {


/// @brief Constructor
///
ResidueArrayAnnealingEvaluator::ResidueArrayAnnealingEvaluator():
	core::pack::interaction_graph::AnnealableGraphBase(),
	weighted_energy_methods_(),
	source_pose_residues_(),
	current_residues_(),
	//current_energy_,
	considered_node_(0),
	considered_state_(0),
	//considered_energy_,
	per_node_rotamer_sets_()
{
	current_energy_ = 0.0;
	considered_energy_ = 0.0;
}

///@brief Destructor
//.
ResidueArrayAnnealingEvaluator::~ResidueArrayAnnealingEvaluator()
{}

///@brief Copy constructor
///
ResidueArrayAnnealingEvaluator::ResidueArrayAnnealingEvaluator( ResidueArrayAnnealingEvaluator const &src ) :
	core::pack::interaction_graph::AnnealableGraphBase( src ),
	weighted_energy_methods_( src.weighted_energy_methods_ ),
	source_pose_residues_( src.source_pose_residues_ ),
	current_residues_( src.current_residues_ ),
	//current_energy_,
	considered_node_( src.considered_node_ ),
	considered_state_( src.considered_state_ ),
	//considered_energy_,
	per_node_rotamer_sets_( src.per_node_rotamer_sets_ )
{
	current_energy_=src.current_energy_;
	considered_energy_=src.considered_energy_;
}

void ResidueArrayAnnealingEvaluator::initialize(
	core::scoring::ScoreFunction const & score_function,
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets & rotamer_sets,
	graph::GraphCOP )
{
	runtime_assert_string_msg(
		rotamer_sets.total_residue() == pose.total_residue(),
		"rotamer_sets.total_residue() != pose.total_residue()" );

	// Get energy method
	foreach_ ( core::scoring::methods::WholeStructureEnergyOP energy_method, std::make_pair(score_function.ws_methods_begin(), score_function.ws_methods_end()) ) {
		core::scoring::annealing::ResidueArrayAnnealableEnergyOP annealable_energy_method ( utility::pointer::dynamic_pointer_cast< core::scoring::annealing::ResidueArrayAnnealableEnergy >( energy_method ) );

		if ( ! annealable_energy_method ) {
			continue;
		}

		runtime_assert_string_msg(
			energy_method->score_types().size() == 1,
			"annealable_energy_method reported multiple supported score types.");

		core::PackerEnergy score_weight = score_function.get_weight( energy_method->score_types().front() );
		if ( score_weight == 0.0 ) {
			continue;
		}

		weighted_energy_methods_.push_back( std::make_pair( score_weight, annealable_energy_method ));
	}

	// Store rotamer sets for annealing
	for ( core::Size i = 1; i <= rotamer_sets.nmoltenres(); ++i ) {
		per_node_rotamer_sets_.push_back(rotamer_sets.rotamer_set_for_moltenresidue( i ));
	}

	// Setup residue arrays for current and considered steps
	source_pose_residues_.resize( rotamer_sets.total_residue() );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		source_pose_residues_[i] = pose.residue(i).get_self_ptr();
	}

	blanket_assign_state_0();
}

int ResidueArrayAnnealingEvaluator::get_num_nodes() const
{
	return per_node_rotamer_sets_.size();
}

int ResidueArrayAnnealingEvaluator::get_num_states_for_node(int n) const
{
	return per_node_rotamer_sets_.at(n)->num_rotamers();
}

int ResidueArrayAnnealingEvaluator::get_num_total_states() const
{
	core::Size count = 0;

	for ( core::Size i = 1; i <= per_node_rotamer_sets_.size(); ++i ) {
		count += per_node_rotamer_sets_[i]->num_rotamers();
	}

	return count;
}

void ResidueArrayAnnealingEvaluator::prepare_for_simulated_annealing()
{
}

void ResidueArrayAnnealingEvaluator::clear_consideration()
{
	considered_node_ = 0;
	considered_state_ = 0;
	considered_energy_ = current_energy_;
}

void ResidueArrayAnnealingEvaluator::blanket_assign_state_0()
{
	// State zero is not well defined for this model, assigned to the pose's origin residue.
	current_residues_ = source_pose_residues_;
	current_energy_ = calculate_weighted_energy( current_residues_ );

	clear_consideration();
}

bool ResidueArrayAnnealingEvaluator::any_vertex_state_unassigned() const
{
	for ( int node_ind = 1; node_ind <= get_num_nodes(); ++node_ind ) {
		int resid = per_node_rotamer_sets_.at(node_ind)->resid();

		if ( current_residues_.at(resid) == source_pose_residues_.at(resid) ) {
			return true;
		}
	}

	return false;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::set_state_for_node(int node_ind , int new_state)
{
	current_residues_.at(per_node_rotamer_sets_.at(node_ind)->resid()) = per_node_rotamer_sets_.at(node_ind)->rotamer(new_state);
	current_energy_ = calculate_weighted_energy( current_residues_ );

	clear_consideration();

	return current_energy_;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::set_network_state( ObjexxFCL::FArray1_int & state_array )
{
	for ( int node_ind = 1; node_ind <= get_num_nodes(); ++node_ind ) {
		int new_state = state_array(node_ind);
		current_residues_.at(per_node_rotamer_sets_.at(node_ind)->resid()) = per_node_rotamer_sets_.at(node_ind)->rotamer(new_state);
	}

	current_energy_ = calculate_weighted_energy( current_residues_ );

	clear_consideration();

	return current_energy_;
}

void ResidueArrayAnnealingEvaluator::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node)
{

	int node_resid = per_node_rotamer_sets_.at(node_ind)->resid();
	core::conformation::ResidueCOP current_res = current_residues_.at(node_resid);

	current_residues_.at( node_resid ) = per_node_rotamer_sets_.at(node_ind)->rotamer(new_state);

	considered_node_ = node_ind;
	considered_state_ = new_state;
	considered_energy_ = calculate_weighted_energy( current_residues_ );

	current_residues_.at( node_resid ) = current_res;

	delta_energy = considered_energy_ - current_energy_;
	//TODO alexford Need previous energy?
	prev_energy_for_node = 0;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::commit_considered_substitution()
{
	int node_resid = per_node_rotamer_sets_.at(considered_node_)->resid();
	current_residues_.at( node_resid ) = per_node_rotamer_sets_.at(considered_node_)->rotamer(considered_state_);
	current_energy_ = considered_energy_;

	return current_energy_;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::get_energy_current_state_assignment()
{
	return current_energy_;
}

core::Real ResidueArrayAnnealingEvaluator::calculate_weighted_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect ) {
	core::Real result = 0;

	foreach_ ( WeightedMethodPair method, weighted_energy_methods_ ) {
		result += method.first * method.second->calculate_energy( resvect );
	}

	return result;
}

void ResidueArrayAnnealingEvaluator::set_errorfull_deltaE_threshold( core::PackerEnergy )
{
}

}
}
}
