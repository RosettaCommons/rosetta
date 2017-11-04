// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

//Symmetry headers:
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>


#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <utility/graph/Graph.fwd.hh>

#include <basic/Tracer.hh>

#include <boost/foreach.hpp>
#define foreach_ BOOST_FOREACH
#include <utility>

static THREAD_LOCAL basic::Tracer TR("core.pack.interaction_graph.ResidueArrayAnnealingEvaluator");

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
	per_node_rotamer_sets_(),
	is_symmetric_( false ),
	has_mirror_symm_( false ),
	num_indep_nodes_(0),
	dependent_node_map_(),
	dependent_residue_map_()
{
	current_energy_ = 0.0;
	considered_energy_ = 0.0;
}

/// @brief Destructor
//.
ResidueArrayAnnealingEvaluator::~ResidueArrayAnnealingEvaluator()
{}

/// @brief Copy constructor
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
	per_node_rotamer_sets_( src.per_node_rotamer_sets_ ),
	is_symmetric_( src.is_symmetric_ ),
	has_mirror_symm_( src.has_mirror_symm_ ),
	num_indep_nodes_( src.num_indep_nodes_),
	dependent_node_map_( src.dependent_node_map_ ),
	dependent_residue_map_( src.dependent_residue_map_ )
{
	current_energy_=src.current_energy_;
	considered_energy_=src.considered_energy_;
}

void ResidueArrayAnnealingEvaluator::initialize(
	core::scoring::ScoreFunction const & score_function,
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets & rotamer_sets,
	utility::graph::GraphCOP )
{
	runtime_assert_string_msg(
		rotamer_sets.total_residue() == pose.size(),
		"rotamer_sets.total_residue() != pose.size()" );

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

		// Call the setup_residuearrayannealableenergy_for_packing() function for the method, which allows
		// the method to cache data from the Pose (which will be invisible to it during packing).
		annealable_energy_method->setup_residuearrayannealableenergy_for_packing( pose, score_function );

		weighted_energy_methods_.push_back( std::make_pair( score_weight, annealable_energy_method ));
	}

	if ( weighted_energy_methods_.empty() ) return; //All of the expensive setup that follows this point can be skipped if there are no energy methods that will use this evaluator.

	// Store rotamer sets for annealing
	for ( core::Size i = 1; i <= rotamer_sets.nmoltenres(); ++i ) {
		per_node_rotamer_sets_.push_back(rotamer_sets.rotamer_set_for_moltenresidue( i ));
	}
	num_indep_nodes_ = per_node_rotamer_sets_.size(); //Count number of indpendent nodes (prior to symmetry setup).

	// Set up the symmetry information (if this is a symmetric packing job):
	initialize_symmetry_info( pose ); //Also expands the per_node_rotamer_sets_ for the symmetric copies.

	// Setup residue arrays for current and considered steps
	source_pose_residues_.resize( rotamer_sets.total_residue() );
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		source_pose_residues_[i] = pose.residue(i).get_self_ptr();
	}

	blanket_assign_state_0();
}

/// @brief Get the number of nodes.
/// @details This lies a little bit.  It only returns the number of independent nodes.
int ResidueArrayAnnealingEvaluator::get_num_nodes() const
{
	return num_indep_nodes_;
}

/// @brief Get the number of states for a specific node.
/// @param[in] node Index of the node.
int ResidueArrayAnnealingEvaluator::get_num_states_for_node(int n) const
{
	return per_node_rotamer_sets_.at(n)->num_rotamers();
}

/// @brief Get the total number of states for all nodes.
/// @details This lies a little bit.  It only returns the total number of states for the independent nodes.
int ResidueArrayAnnealingEvaluator::get_num_total_states() const
{
	core::Size count = 0;

	for ( core::Size i = 1; i <= num_indep_nodes_; ++i ) { //Only iterate over the independent nodes when counting states.
		count += per_node_rotamer_sets_[i]->num_rotamers();
	}

	return count;
}

void ResidueArrayAnnealingEvaluator::prepare_for_simulated_annealing()
{
}

/// @brief Initialize symmetry information.
/// @details Called by initialize().
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueArrayAnnealingEvaluator::initialize_symmetry_info(
	core::pose::Pose const &pose
) {
	debug_assert( num_indep_nodes_ == per_node_rotamer_sets_.size() );

	core::conformation::ConformationCOP conf( pose.conformation_ptr() );
	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( conf ) ); // Will be NULL if pose is asymmetric.
	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorsymmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const >( conf ) ); // Will be NULL if pose is asymmetric or if there are no mirror operations.
	is_symmetric_ = static_cast<bool>( symmconf );
	has_mirror_symm_ = static_cast<bool>( mirrorsymmconf );

	dependent_node_map_.clear();
	dependent_residue_map_.clear();

	if ( !is_symmetric_ ) return;
	core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
	for ( core::Size inode=1, inodemax=num_indep_nodes_; inode<=inodemax; ++inode ) { //Loop through the nodes that are already present in the per_node_rotamer_sets_ vector.
		core::Size const curres( per_node_rotamer_sets_.at(inode)->resid() ); //The independent residue index.
		debug_assert( symminfo->bb_is_independent(curres) ); //Should be true.
		utility::vector1 < core::Size > child_node_indices;
		utility::vector1 < core::Size > child_residue_indices;
		for ( core::Size ir=1, irmax=symminfo->num_total_residues_without_pseudo(); ir<=irmax; ++ir ) { //Loop through all pose residues
			if ( !symminfo->bb_is_independent(ir) ) {
				core::Size const parent( symminfo->bb_follows(ir) );
				if ( parent == curres ) {
					core::pack::rotamer_set::symmetry::SymmetricRotamerSet_COP symrotset( utility::pointer::dynamic_pointer_cast< core::pack::rotamer_set::symmetry::SymmetricRotamerSet_ const >( per_node_rotamer_sets_.at(inode) ) );
					runtime_assert( symrotset ); //Should be a SymetricRotamerSet_.
					per_node_rotamer_sets_.push_back( symrotset->orient_rotamer_set_to_symmetric_partner( pose, pose.residue(curres), ir, *symrotset, has_mirror_symm_ ) /*Creates a copy and returns OP.  Should be mirror-compatible.*/ );
					//update_residue_index( per_node_rotamer_sets_[per_node_rotmaer_sets_.size(), ir);
					child_node_indices.push_back( per_node_rotamer_sets_.size() );
					child_residue_indices.push_back( ir );
				}
			}
		}
		dependent_node_map_[ static_cast<int>(inode) ] = child_node_indices;
		dependent_residue_map_[ static_cast<int>(inode) ] = child_residue_indices;
	}

	//Debug output:
	if ( TR.Debug.visible() ) {
		TR.Debug << "Completed ResidueArrayAnnealingEvaluator::initialize_symmetry_info()." << std::endl;
		TR.Debug << "Node\tDepNodes\tDepResidues" << std::endl;
		for ( core::Size i=1; i<=num_indep_nodes_; ++i ) {
			TR.Debug << i << "\t";
			for ( core::Size j=1, jmax=dependent_node_map_.at(i).size(); j<=jmax; ++j ) {
				TR.Debug << dependent_node_map_.at(i)[j];
				if ( j<jmax ) TR.Debug << ",";
			}
			TR.Debug << "\t";
			for ( core::Size j=1, jmax=dependent_residue_map_.at(i).size(); j<=jmax; ++j ) {
				TR.Debug << dependent_residue_map_.at(i)[j];
				if ( j<jmax ) TR.Debug << ",";
			}
			TR.Debug << std::endl;
		}
		TR.flush();
	}

}

/// @brief Sets the current consideration.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void ResidueArrayAnnealingEvaluator::set_consideration (
	int const node_ind,
	int const node_resid,
	int const new_state,
	utility::vector1 < std::pair < int, core::conformation::ResidueCOP > > & unset_information
) {
	using namespace core::conformation;

	unset_information.push_back( std::pair< int, ResidueCOP>( node_resid, current_residues_.at(node_resid) ) );

	current_residues_.at( node_resid ) = per_node_rotamer_sets_.at(node_ind)->rotamer(new_state);

	//If this is the symmetric case, we also need to update the symmetric nodes:
	if ( is_symmetric_ ) {
		utility::vector1 < core::Size > const & dependent_nodes( dependent_node_map_.at( node_ind ) );
		utility::vector1 < core::Size > const & dependent_residues( dependent_residue_map_.at( node_ind ) );
		for ( core::Size i=1, imax=dependent_nodes.size(); i<=imax; ++i ) {
			unset_information.push_back( std::pair< int, ResidueCOP>( static_cast<int>( dependent_residues[i] ), current_residues_.at( dependent_residues[i] ) ) );
			current_residues_.at( dependent_residues[i] ) = per_node_rotamer_sets_.at( dependent_nodes[i] )->rotamer(new_state);
		}
	}
}

/// @brief Clears the current consideration.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void ResidueArrayAnnealingEvaluator::unset_consideration(
	utility::vector1< std::pair < int, core::conformation::ResidueCOP> > const & unset_info
) {
	for ( core::Size i=1, imax=unset_info.size(); i<=imax; ++i ) {
		current_residues_.at( unset_info[i].first ) = unset_info[i].second;
	}
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

	//If this is the symmetric case, we also need to update the symmetric nodes:
	if ( is_symmetric_ ) {
		utility::vector1 < core::Size > const & dependent_nodes( dependent_node_map_.at( node_ind ) );
		utility::vector1 < core::Size > const & dependent_residues( dependent_residue_map_.at( node_ind ) );
		for ( core::Size i=1, imax=dependent_nodes.size(); i<=imax; ++i ) {
			current_residues_.at( dependent_residues[i] ) = per_node_rotamer_sets_.at( dependent_nodes[i] )->rotamer(new_state);
		}
	}

	current_energy_ = calculate_weighted_energy( current_residues_ );

	clear_consideration();

	return current_energy_;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::set_network_state( ObjexxFCL::FArray1_int & state_array )
{
	for ( int node_ind = 1; node_ind <= get_num_nodes(); ++node_ind ) {
		int new_state = state_array(node_ind);
		current_residues_.at(per_node_rotamer_sets_.at(node_ind)->resid()) = per_node_rotamer_sets_.at(node_ind)->rotamer(new_state);

		//If this is the symmetric case, we also need to update the symmetric nodes:
		if ( is_symmetric_ ) {
			utility::vector1 < core::Size > const & dependent_nodes( dependent_node_map_.at( node_ind ) );
			utility::vector1 < core::Size > const & dependent_residues( dependent_residue_map_.at( node_ind ) );
			for ( core::Size i=1, imax=dependent_nodes.size(); i<=imax; ++i ) {
				current_residues_.at( dependent_residues[i] ) = per_node_rotamer_sets_.at( dependent_nodes[i] )->rotamer(new_state);
			}
		}
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

	utility::vector1 < std::pair < int, core::conformation::ResidueCOP > > unset_information;
	set_consideration( node_ind, node_resid, new_state, unset_information );

	considered_node_ = node_ind;
	considered_state_ = new_state;
	considered_energy_ = calculate_weighted_energy( current_residues_, considered_node_ );

	unset_consideration( unset_information );

	delta_energy = considered_energy_ - current_energy_;
	//TODO alexford Need previous energy?
	prev_energy_for_node = 0;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::commit_considered_substitution()
{
	for ( WeightedMethodPair const &method : weighted_energy_methods_ ) {
		method.second->commit_considered_substitution();
	}

	int node_resid = per_node_rotamer_sets_.at(considered_node_)->resid();

	utility::vector1<std::pair < int, core::conformation::ResidueCOP > > unset_information_dummy;
	set_consideration( considered_node_, node_resid, considered_state_, unset_information_dummy);
	current_energy_ = considered_energy_;

	return current_energy_;
}

core::PackerEnergy ResidueArrayAnnealingEvaluator::get_energy_current_state_assignment()
{
	return current_energy_;
}

core::Real ResidueArrayAnnealingEvaluator::calculate_weighted_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, int const substitution_position ) {
	core::Real result = 0;

	for ( WeightedMethodPair const & method : weighted_energy_methods_ ) {
		result += method.first * method.second->calculate_energy( resvect, substitution_position );
	}

	return result;
}

void ResidueArrayAnnealingEvaluator::set_errorfull_deltaE_threshold( core::PackerEnergy )
{
}

}
}
}
