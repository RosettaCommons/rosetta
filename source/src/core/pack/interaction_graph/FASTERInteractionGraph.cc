// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

// Utility headers
#include <utility/string_util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>

// STL Headers
#include <list>
#include <algorithm>
#include <iostream>
#include <utility/assert.hh>

// Numeric headers
#include <numeric/random/random.hh>


#include <utility/vector1.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {


/// @brief main constructor, no default or copy constructors
///
/// allocates one-body energy array and initializes it to zero.
///
FASTERNode::FASTERNode(
	InteractionGraphBase * owner,
	int node_id,
	int num_states
) :
	PrecomputedPairEnergiesNode( owner, node_id, num_states ),
	one_body_energies_(num_states + 1, core::PackerEnergy( 0.0 )),
	curr_state_total_energy_( 0.0 ),
	alternate_state_is_being_considered_( false ),
	have_prepared_once_for_FASTER_( false ),
	in_FASTER_mode_( false ),
	perturbed_( false ),
	have_relaxed_since_neighbors_perturbation_( true ),
	have_contributed_deltaE_following_perturbation_( true ),
	relaxed_state_( 0 )
{}

/// @brief destructor
///
/// not responsible for any dynamically allocated memory, so node does nothing
/// it's member variables, of course, are implicitly destructed
FASTERNode::~FASTERNode() = default;

/// @brief prints a description of the node and all of it's one-body energies
void FASTERNode::print() const
{
	std::cout << "NODE: " << get_node_index() << " with " <<
		get_num_states() << " states" << std::endl;
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		std::cout << "(" << ii << ", ";
		std::cout << one_body_energies_[ ii ] << ") ";
		if ( ii % 3 == 0 ) std::cout << std::endl;
	}
	std::cout << std::endl  << "-----------------" << std::endl;
}


/// @brief update energy to the one-body energy for state
///
///
/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be set.
void FASTERNode::update_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
	return;
}

/// @brief set all the one-body energies for this node
///
/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void FASTERNode::update_one_body_energies( FArray1< core::PackerEnergy > & energies )
{
	debug_assert( energies.size() == (unsigned int) get_num_states() );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] = energies( ii );
	}
	return;
}

/// @brief adds energy to the one-body energy for state state
///
/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be added.
void FASTERNode::add_to_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] += energy;
	return;
}

/// @brief adds all the energies in energies to the one-body energies for this node
///
/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void FASTERNode::add_to_one_body_energies( FArray1< core::PackerEnergy > & energies )
{
	debug_assert( energies.size() == (unsigned int) get_num_states() );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] += energies( ii );
	}
	return;
}

/// @brief sets all of the one-body energies for this node to zero
void FASTERNode::zero_one_body_energies()
{
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] = core::PackerEnergy( 0.0 );
	}
}

/// @brief returns the one body energy for a state
///
/// @param state - [in]
core::PackerEnergy FASTERNode::get_one_body_energy( int state )
{
	return one_body_energies_[ state ];
}

/// @brief prepares node for simulated annealing
///
/// updates internal edge vector + other vectorized edge information
void FASTERNode::prepare_for_simulated_annealing()
{
	if ( ! get_edge_vector_up_to_date() ) update_internal_vectors();
	in_FASTER_mode_ = false;
	return;
}

void
FASTERNode::prepare_for_FASTER()
{
	if ( ! get_edge_vector_up_to_date() ) update_internal_vectors();
	in_FASTER_mode_ = true;

	//allocate arrays if necessary;
	if ( ! have_prepared_once_for_FASTER_ ) {
		state_energies_in_current_state_assignment_.resize( get_num_states() + 1 );
		state_energies_in_current_context_.resize( get_num_states() + 1 );
		neighbor_relaxed_in_sBR_.resize( get_num_incident_edges() + 1 );
		perturbed_two_body_energies_.resize( get_num_incident_edges() + 1 );
		have_prepared_once_for_FASTER_ = true;
	}

	//get all rotamer's total energies for current state assignment
	get_total_energy_in_curr_state_assignment_for_all_states();
}

/// @brief assigns node's state to it's zero, or "unassigned" state.
///
/// zeros the edge-energy array, informs neighbors that it's in its unassigned
/// state
void FASTERNode::assign_zero_state()
{

	//std::cout << "assign_state: node -  " << get_node_index() <<
	// " new state " << 0 << "...";

	current_state_ = 0;
	relaxed_state_ = 0;
	alternate_state_ = 0;
	alternate_state_is_being_considered_ = false;

	curr_state_one_body_energy_ = core::PackerEnergy( 0.0 );
	//fills from [1] to end
	auto position1 = curr_state_two_body_energies_.begin();
	++position1;
	std::fill( position1, curr_state_two_body_energies_.end(), core::PackerEnergy( 0.0 ));
	curr_state_total_energy_ = core::PackerEnergy( 0.0 );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_faster_edge(ii)->acknowledge_state_zeroed( get_node_index() );
	}

	return;
}


/// @brief assigns node a new_state
///
/// node updates its curr_state one and two body energies
///
/// @param new_state - [in] - the new state the node should be assigned
void FASTERNode::assign_state(int new_state)
{
	debug_assert( new_state >= 0 && new_state <= get_num_states());

	if ( new_state == 0 ) {
		assign_zero_state();
	} else {
		//std::cout << "assign_state: node -  " << get_node_index() <<
		// " new state " << new_state << "...";
		current_state_ = new_state;
		relaxed_state_ = new_state;
		curr_state_one_body_energy_ = one_body_energies_[ current_state_ ];
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_faster_edge(ii)->acknowledge_state_change(
				get_node_index(),
				current_state_,
				curr_state_two_body_energies_[ii]);

			curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
		}
		//std::cout<< "..done" << std::endl;
	}
	return;
}

/// @brief returns the state the node is currently assigned
int FASTERNode::get_current_state() const
{
	return current_state_;
}

/// @brief returns the one body energy for the state the node is currently assigned
core::PackerEnergy FASTERNode::get_one_body_energy_current_state() const
{ return curr_state_one_body_energy_; }

/// @brief tells the node that it should change its state to the last state it was
/// asked to consider (from a call to project_deltaE_for_substitution)
///
/// updates edge energy vector, iterates across neighbors having them update
/// their edge energies.  Bookkeeping recaptures performance lost by
/// leaving energy2b structure
void FASTERNode::commit_considered_substitution()
{
	debug_assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	relaxed_state_ = current_state_;
	curr_state_one_body_energy_ = alternate_state_one_body_energy_;
	curr_state_total_energy_ = alternate_state_total_energy_;

	//copies from [1] to end
	auto alt_position1 = alternate_state_two_body_energies_.begin();
	++alt_position1;
	auto curr_position1 = curr_state_two_body_energies_.begin();
	++curr_position1;

	std::copy( alt_position1, alternate_state_two_body_energies_.end(), curr_position1 );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_faster_edge(ii)->acknowledge_substitution(
			get_node_index(),
			alternate_state_two_body_energies_[ii],
			current_state_
		);
	}

	alternate_state_is_being_considered_ = false;
	return;
}


/// @brief outputs to standard error the bookkeeping energies for the node in its
/// current state assignment
void FASTERNode::print_internal_energies() const
{
	std::cout << "curr_state " << current_state_ << " ";
	std::cout << "curr_state_one_body_energy_ ";
	std::cout << curr_state_one_body_energy_ << " ";
	std::cout << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		std::cout << "(" << get_index_of_adjacent_node( ii ) << ": " <<
			curr_state_two_body_energies_[ ii ] << ") ";
	}
	std::cout << std::endl;
}

/// @brief removes numerical drift long stretches of efficient bookkeeping
/// produces
void FASTERNode::update_internal_energy_sums()
{
	debug_assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = core::PackerEnergy( 0.0 );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_total_energy_ += get_incident_faster_edge(ii)->get_current_two_body_energy();
	}
	curr_state_total_energy_ += curr_state_one_body_energy_;
	return;
}

/// @brief If FASTERNode is the most-derived class being used, then this function
/// will be called and will return the amount of memory statically allocated by
/// a single FASTERNode.
unsigned int
FASTERNode::count_static_memory() const
{
	return sizeof( FASTERNode );
}

/// @brief Called either by the IGBase if the FASTERNode is the most-derived class, or
/// called recursively by a derived class.  Called to account for the dynamically allocated
/// memory that this node uses.
unsigned int
FASTERNode::count_dynamic_memory() const
{
	unsigned int dynamic_memory = 0;
	dynamic_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );
	dynamic_memory += neighbors_curr_state_.size() * sizeof( int );
	dynamic_memory += edge_matrix_ptrs_.size() * sizeof ( FArray2< core::PackerEnergy > );
	dynamic_memory += curr_state_two_body_energies_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += alternate_state_two_body_energies_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += state_energies_in_current_state_assignment_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += state_energies_in_current_context_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += neighbor_relaxed_in_sBR_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += perturbed_two_body_energies_.size() * sizeof ( core::PackerEnergy );

	dynamic_memory += NodeBase::count_dynamic_memory();
	return dynamic_memory;
}

/// @brief updates bookkeeping arrays that correspond to edge-list.
///
/// calls base class update_edge_vector function, and then proceeds to create
/// appropriate bookkeeping arrays used in simulated annealing
void FASTERNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	//aa_offsets_for_this_lookup_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_.resize( get_num_incident_edges() + 1);

	edge_matrix_ptrs_.clear();
	edge_matrix_ptrs_.reserve( get_num_incident_edges() + 1);
	edge_matrix_ptrs_.emplace_back( ); //occupy the 0th position

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		edge_matrix_ptrs_.push_back( get_incident_faster_edge(ii)->get_edge_table_ptr() );
	}

	curr_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	alternate_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	return;
}


void
FASTERNode::get_total_energy_in_curr_state_assignment_for_all_states()
{
	std::copy(one_body_energies_.begin(),
		one_body_energies_.end(),
		state_energies_in_current_state_assignment_.begin()
	);
	int const local_num_states = get_num_states();

	for ( int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes(); ++ii ) {

		if ( neighbors_curr_state_[ ii ] == 0 ) continue;
		// Lower neighbor;
		// since the edge table is allocated( node2_nstates, node1_nstates )
		// and we're "node 2", we can walk across the row very efficiently

		FArray2A< core::PackerEnergy > const & edge_table( edge_matrix_ptrs_[ ii ] );
		int li_curr = edge_table.index( 1, neighbors_curr_state_[ ii ] );
		for ( int jj = 1; jj <= local_num_states; ++jj, ++li_curr ) {
			state_energies_in_current_state_assignment_[ jj ] += edge_table[ li_curr ];
		}
	}
	for ( int ii = get_num_edges_to_smaller_indexed_nodes() + 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( neighbors_curr_state_[ ii ] == 0 ) continue;

		// Upper neighbor;
		// Since the edge table is allocated( node2_nstates, node1_nstates )
		// and we're "node 1", we have to take strides across the table with the
		// step size equal to the number of states on the other node.

		FArray2A< core::PackerEnergy > const & edge_table( edge_matrix_ptrs_[ ii ] );
		int const stride = edge_table.size1();
		int li_curr = edge_table.index( neighbors_curr_state_[ ii ], 1 );
		for ( int jj = 1; jj <= local_num_states; ++jj, li_curr += stride ) {
			state_energies_in_current_state_assignment_[ jj ] += edge_table[ li_curr ];
		}
	}
}

void FASTERNode::partial_assign_state_with_lowest_one_body_energy()
{
	core::PackerEnergy best_one_body_energy( 12345 );
	int state_with_best_one_body_energy( 0 );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		if ( ii == 1 || one_body_energies_[ ii ] < best_one_body_energy ) {
			best_one_body_energy = one_body_energies_[ ii ];
			state_with_best_one_body_energy = ii;
		}
	}
	partial_assign_state( state_with_best_one_body_energy );
}

void FASTERNode::partial_assign_relaxed_state( Real prob )
{
	bool accept = prob < 1.0 ? (numeric::random::rg().uniform() < prob) : true;

	if ( accept ) {
		partial_assign_state( relaxed_state_ );
	} else {
		partial_assign_state( current_state_ );
	}
}

void FASTERNode::partial_assign_state( int new_state )
{
	current_state_ = new_state;
	relaxed_state_ = new_state;
	curr_state_total_energy_ = curr_state_one_body_energy_ =
		one_body_energies_[ current_state_ ];

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_faster_edge( ii )->acknowledge_partial_state_assignment(
			get_node_index(),
			current_state_
		);
	}
}

void FASTERNode::acknowledge_neighbors_partial_state_assignment(
	int which_neighbor,
	int neighbors_new_state
)
{
	neighbors_curr_state_[ which_neighbor ] = neighbors_new_state;
}

void FASTERNode::complete_partial_state_assignment()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_two_body_energies_[ ii ] =
			get_incident_faster_edge( ii )->get_curr_state_energy_following_partial_state_assignment();
		curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	}

	if ( in_FASTER_mode_ ) {
		get_total_energy_in_curr_state_assignment_for_all_states();
	}
}

void FASTERNode::acknowledge_neighbors_perturbed_state(
	int which_neighbor,
	int const neighbors_perturbed_state
)
{
	int const local_num_states = get_num_states();

	if ( perturbed_ ) return;

	have_relaxed_since_neighbors_perturbation_ = false;
	FArray2A< core::PackerEnergy > const & edge_table( edge_matrix_ptrs_[ which_neighbor ] );
	/*int li_curr = edge_table.index( 1, neighbors_curr_state_[ which_neighbor ] );
	int li_pert = edge_table.index( 1, perturbed_state );

	/// This loop is the most time consuming step of FASTER -- it is made as minimal as possible.
	for ( int ii = 1; ii <= local_num_states; ++ii, ++li_curr, ++li_pert ) {
	state_energies_in_current_context_[ ii ] += edge_table[ li_pert ] - edge_table[ li_curr ];
	}*/

	int const neighbors_current_state = neighbors_curr_state_[ which_neighbor ];
	if ( which_neighbor <= get_num_edges_to_smaller_indexed_nodes() ) {
		/*for ( int ii = 1; ii <= local_num_states; ++ii ) {
		state_energies_in_current_context_[ ii ] +=
		get_incident_faster_edge( which_neighbor )->get_two_body_energy(
		neighbors_perturbed_state, ii )
		-
		get_incident_faster_edge( which_neighbor )->get_two_body_energy(
		neighbors_curr_state, ii );

		}*/

		// Lower neighbor;
		// since the edge table is allocated( node2_nstates, node1_nstates )
		// and we're "node 2", we can walk across the row very efficiently
		int li_curr = edge_table.index( 1, neighbors_current_state );
		int li_pert = edge_table.index( 1, neighbors_perturbed_state );

		/// This loop is the most time consuming step of FASTER -- it is made as minimal as possible.
		for ( int ii = 1; ii <= local_num_states; ++ii, ++li_curr, ++li_pert ) {
			state_energies_in_current_context_[ ii ] += edge_table[ li_pert ] - edge_table[ li_curr ];
		}
	} else {
		/*for ( int ii = 1; ii <= local_num_states; ++ii ) {
		state_energies_in_current_context_[ ii ] +=
		get_incident_faster_edge( which_neighbor )->get_two_body_energy(
		ii, neighbors_perturbed_state );
		-
		get_incident_faster_edge( which_neighbor )->get_two_body_energy(
		ii, neighbors_current_state);
		}*/

		// Upper neighbor;
		// Since the edge table is allocated( node2_nstates, node1_nstates )
		// and we're "node 1", we have to take strides across the table with the
		// step size equal to the number of states on the other node.
		int const stride = edge_table.size1();
		int li_curr = edge_table.index( neighbors_current_state, 1 );
		int li_pert = edge_table.index( neighbors_perturbed_state, 1 );
		for ( int ii = 1; ii <= local_num_states; ++ii, li_curr += stride, li_pert += stride ) {
			state_energies_in_current_context_[ ii ] += edge_table[ li_pert ] - edge_table[ li_curr ];
		}

	}

}

void
FASTERNode::prepare_for_perturbation()
{
	perturbed_ = true;
	have_contributed_deltaE_following_perturbation_ = false;
}

void
FASTERNode::set_no_longer_perturbed()
{
	perturbed_ = false;
}

void
FASTERNode::set_perturbed_state( int perturbed_state )
{
	relaxed_state_ = perturbed_state;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		perturbed_two_body_energies_[ ii ] = get_incident_faster_edge( ii )->
			get_energy_for_perturbed_state( get_node_index(), relaxed_state_ );
		get_incident_faster_edge( ii )->acknowledge_participation_in_perturbation();
	}
}

void FASTERNode::relax_neighbors()
{
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( neighbor_relaxed_in_sBR_[ ii ] ) {
			get_incident_faster_edge( ii )->acknowledge_perturbed_state( get_node_index(), relaxed_state_ );
			get_adjacent_faster_node( ii )->relax_after_neighbors_perturbation();
		}
	}

}

void FASTERNode::relax()
{
	int const local_num_states = get_num_states();

	core::PackerEnergy best_relaxed_energy = state_energies_in_current_context_[ 1 ];
	relaxed_state_ = 1;
	for ( int ii = 2; ii <= local_num_states; ++ii ) {
		if ( state_energies_in_current_context_[ ii ] < best_relaxed_energy ) {
			best_relaxed_energy = state_energies_in_current_context_[ ii ];
			relaxed_state_ = ii;
		}
	}

}

void FASTERNode::reset_relaxed_for_neighbors()
{
	reset_relaxed();
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_adjacent_faster_node( ii )->reset_relaxed();
	}
}


void FASTERNode::reset_relaxed()
{
	relaxed_state_ = current_state_;
}


void FASTERNode::tell_neighbors_to_prep_for_relaxation()
{
	int const num_to_relax = 10;


	core::PackerEnergy top10[ num_to_relax ];
	int which10[ num_to_relax ];
	for ( int ii = 0; ii < num_to_relax; ++ii ) {
		top10[ ii ] = 0;
		which10[ ii ] = 0;
	}

	//decide which neighbors to relax
	if ( get_num_incident_edges() > num_to_relax ) {
		int whichBottom;
		core::PackerEnergy bottom_of_top10;
		for ( int ii = 1; ii <= num_to_relax; ++ii ) {
			top10[ ii - 1 ] = std::abs(perturbed_two_body_energies_[ ii ]);
			which10[ ii - 1] = ii;
			if ( ii == 1 || bottom_of_top10 > top10[ ii-1 ] ) {
				bottom_of_top10 = top10[ ii-1 ];
				whichBottom = ii-1;
			}
		}

		for ( int ii = 11; ii <= get_num_incident_edges(); ++ii ) {
			core::PackerEnergy iiabs = std::abs( perturbed_two_body_energies_[ ii ] );
			if ( iiabs > bottom_of_top10 ) {
				top10[ whichBottom ] = iiabs;
				which10[ whichBottom ] = ii;
				for ( int jj = 0; jj < num_to_relax; ++jj ) {
					if ( jj == 0 || bottom_of_top10 > top10[ jj ] ) {
						whichBottom = jj;
						bottom_of_top10 = top10[ jj ];
					}
				}
			}
		}
		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			neighbor_relaxed_in_sBR_[ ii ] = false;
		}
		//std::cout << "neighbors: ";
		for ( int ii : which10 ) {
			//std::cout << which10[ ii ] << " " ;
			neighbor_relaxed_in_sBR_[ ii ] = true;
		}
		//std::cout << std::endl;
	} else {
		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			neighbor_relaxed_in_sBR_[ ii ] = true;
		}
	}


	//std::cout << "peturbed: " << get_node_index() << " relaxing: ";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( neighbor_relaxed_in_sBR_[ ii ] ) {
			get_adjacent_faster_node( ii )->prep_for_neighbors_perturbation();
			//std::cout << get_index_of_adjacent_node( ii ) << " ";
		}
	}
	//std::cout << std::endl;
}

void FASTERNode::prep_for_neighbors_perturbation()
{
	if ( perturbed_ ) return;

	if ( have_contributed_deltaE_following_perturbation_ ) {
		have_contributed_deltaE_following_perturbation_ = false;

		std::copy( state_energies_in_current_state_assignment_.begin(),
			state_energies_in_current_state_assignment_.end(),
			state_energies_in_current_context_.begin()
		);
		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_faster_edge( ii )->acknowledge_participation_in_perturbation();
		}
	}

}

void FASTERNode::relax_after_neighbors_perturbation()
{
	if ( have_relaxed_since_neighbors_perturbation_ || perturbed_ ) return;
	have_relaxed_since_neighbors_perturbation_ = true;

	int const local_num_states = get_num_states();
	core::PackerEnergy best_relaxed_energy = state_energies_in_current_context_[ 1 ];
	relaxed_state_ = 1;
	for ( int ii = 2; ii <= local_num_states; ++ii ) {
		if ( state_energies_in_current_context_[ ii ] < best_relaxed_energy ) {
			best_relaxed_energy = state_energies_in_current_context_[ ii ];
			relaxed_state_ = ii;
		}
	}
}

core::PackerEnergy
FASTERNode::get_deltaE_for_relaxed_state_following_perturbation()
{
	if ( have_contributed_deltaE_following_perturbation_ ) return core::PackerEnergy( 0.0 );

	have_contributed_deltaE_following_perturbation_ = true;

	auto deltaE = core::PackerEnergy( 0.0 );
	if ( perturbed_ ) {
		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			deltaE += get_incident_faster_edge( ii )->get_deltaE_for_perturbation();

			if ( neighbor_relaxed_in_sBR_[ ii ] ) {
				deltaE += get_incident_faster_edge( ii )->
					get_deltaE_for_neighbor_following_perturbation( get_node_index() );
				neighbor_relaxed_in_sBR_[ ii ] = false;
			}
		}
	} else {
		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			deltaE += get_incident_faster_edge( ii )->get_deltaE_for_perturbation();
		}
	}
	deltaE += ( one_body_energies_[ relaxed_state_ ] - curr_state_one_body_energy_ );
	return deltaE;
}


core::PackerEnergy FASTERNode::get_one_body_energy_for_relaxed_state() const
{
	return one_body_energies_[ relaxed_state_ ];
}

int
FASTERNode::get_relaxed_state() const
{
	return relaxed_state_;
}


int FASTERNode::get_random_neighbor()
{
	if ( get_num_incident_edges() == 0 ) return 0;

	int ran_neighbor = ((int) (numeric::random::rg().uniform() * get_num_incident_edges() )) + 1;
	return get_index_of_adjacent_node( ran_neighbor );
}


/// @brief main constructor - no default nor copy constructors provided
///
/// @param owner - [in] - pointer to the graph that created this node
/// @param first_node_ind - [in] - the index of the smaller-indexed node
/// @param second_node_ind - [in] - the index of the larger-indexed node
FASTEREdge::FASTEREdge(
	InteractionGraphBase * owner,
	int first_node_ind,
	int second_node_ind
) :
	parent( owner, first_node_ind, second_node_ind),
	two_body_energies_(
	get_faster_node(1)->get_num_states(),
	get_faster_node(0)->get_num_states(),
	core::PackerEnergy( 0.0 )
	),
	curr_state_energy_( core::PackerEnergy( 0.0 )),
	energies_updated_since_last_prep_for_simA_( true ),
	have_contributed_deltaE_following_perturbation_( false ),
	partial_state_assignment_( false )
{
}

/// @brief destructor.  All dynamically allocated memory is managed by the objects contained
/// inside the FASTEREdge, so there is no work to be (explicitly) done.
FASTEREdge::~FASTEREdge() = default;

/// @brief adds the input energy to the two body energy for state1 on the node with the
/// smaller index and state2 on the node with the larger index.
void FASTEREdge::add_to_two_body_energy
(
	int const state1,
	int const state2,
	core::PackerEnergy const energy
)
{
	two_body_energies_(state2, state1) += edge_weight() * energy;
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief Adds all the energies stored in the oversized_res_res_energy array to the
/// two body energy table for those states whose amion acid types were
/// previoudsly declared to be amino-acid neighbors.  The res-res array
/// should have the dimension (node1->get_num_states() x node2->get_num_states());
///
/// @param res_res_energy_array - [in] - an array containing the state pair energies
void FASTEREdge::add_to_two_body_energies
(
	FArray2< core::PackerEnergy > const & res_res_energy_array
)
{
	debug_assert( res_res_energy_array.size1() == two_body_energies_.size1() );
	debug_assert( res_res_energy_array.size2() == two_body_energies_.size2() );
	for ( Size ii = 1, iie = two_body_energies_.size2(); ii <= iie; ++ii ) {
		for ( Size jj = 1, jje = two_body_energies_.size1(); jj <= jje; ++jj ) {
			two_body_energies_( jj, ii ) += edge_weight() * res_res_energy_array( jj, ii );
		}
	}
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief Sets the two-body energy for a pair of states.  That is, it overwrites
/// whatever two-body energy there was previously for that state pair with
/// a new energy.  Ignores non-neighboring state pairs.
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
/// @param energy - [in] - the energy which replaces the old two-body energy
void FASTEREdge::set_two_body_energy
(
	int const state1,
	int const state2,
	core::PackerEnergy const energy
)
{
	two_body_energies_( state2, state1 ) = edge_weight() * energy;
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief Sets the two-body energy for a pair of states.  That is, it overwrites
/// whatever two-body energy there was previously for that state pair with
/// a new energy.  Ignores non-neighboring state pairs.
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
/// @param energy - [in] - the energy which replaces the old two-body energy
void FASTEREdge::clear_two_body_energy
(
	int const state1,
	int const state2
)
{
	two_body_energies_(state2,state1) = core::PackerEnergy( 0.0 );
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief returns the two body energy for a pair of states: 0 if those states are
/// not neighbors
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
core::PackerEnergy FASTEREdge::get_two_body_energy( int const state1, int const state2 ) const
{
	return two_body_energies_(state2, state1);
}

/// @brief If all of the energies for an edge have been added in, then declare the edge energies
/// final.  This may mean that the edge deletes itself.
void FASTEREdge::declare_energies_final()
{
	prepare_for_simulated_annealing();
}

/// @brief looks at all pair energies, and if they are all 0, deletes itself
void FASTEREdge::prepare_for_simulated_annealing()
{
	if ( ! energies_updated_since_last_prep_for_simA_ ) return;


	energies_updated_since_last_prep_for_simA_ = false;

	bool any_non_zero = false;
	Size const num_energies = two_body_energies_.size();
	for ( Size ii = 0; ii < num_energies; ++ii ) {
		if ( two_body_energies_[ ii ] != core::PackerEnergy( 0.0 ) ) { any_non_zero = true; break;}
	}

	if ( ! any_non_zero ) delete this;
}

void FASTEREdge::prepare_for_FASTER()
{
	prepare_for_simulated_annealing();
}


/// @brief returns the two body energy corresponding to the current states assigned to
/// the nodes this edge is incident upon.
core::PackerEnergy FASTEREdge::get_current_two_body_energy()
{
	return curr_state_energy_;

}

void FASTEREdge::acknowledge_partial_state_assignment(
	int node,
	int new_state
)
{
	int const other = node == get_node_index(0) ? 1 : 0;

	partial_state_assignment_ = true;
	get_faster_node( other )->acknowledge_neighbors_partial_state_assignment(
		get_edges_position_in_nodes_edge_vector( other ),
		new_state
	);
}

/// @brief updates bookkeeping information when one of the two nodes changes its state
///
/// @param node_ind - [in] - the index of the node that changed its state
/// @param node_state - [in] - the index of the new state it assumed
/// @param new_energy - [out] - the two body energy produced  by the new state and
///  the current state on the other node
void
FASTEREdge::acknowledge_state_change
(
	int node_ind,
	int new_state,
	core::PackerEnergy & new_energy
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	int nodes_curr_states[2];

	nodes_curr_states[ node_substituted ] = new_state;

	nodes_curr_states[ node_not_substituted ] =
		get_faster_node( node_not_substituted )->get_current_state();

	bool one_node_in_zero_state = ( nodes_curr_states[0] == 0 || nodes_curr_states[1] == 0 );

	if (  one_node_in_zero_state ) {
		curr_state_energy_ = core::PackerEnergy( 0.0 );
	} else {
		curr_state_energy_ = two_body_energies_( nodes_curr_states[ 1 ], nodes_curr_states[ 0 ] );
	}
	new_energy = curr_state_energy_;

	get_faster_node( node_not_substituted )->acknowledge_neighbors_state_substitution (
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		new_state
	);

	return;
}

/// @brief updates bookkeeping information when one of the two nodes enters its
/// "unassigned" state.
///
/// @param node_ind - [in] - the index of the node that has just entered its 0 state
void FASTEREdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted = ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = core::PackerEnergy( 0.0 );

	get_faster_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		0
	);
}


/// @brief Returns a reference to the first element in the dense two-body energy
/// table.  Used to create a proxy array on the nodes for cache efficiency.

FArray2A< core::PackerEnergy >
FASTEREdge::get_edge_table_ptr() {
	return FArray2A< core::PackerEnergy >( two_body_energies_( 1, 1 ),
		get_num_states_for_node( 1 ),
		get_num_states_for_node( 0 ) );
}


/// @brief returns the memory usage of the two body energy table for this edge
int FASTEREdge::get_two_body_table_size() const
{
	return two_body_energies_.size();
}

/// @brief returns sizeof FASTEREdge if this is the most-derived instance of the class
unsigned int
FASTEREdge::count_static_memory() const
{
	return sizeof( FASTEREdge );
}

/// @brief returns the amount of memory dynamically allocated by the edge and
/// recurses on its parent (in this case, the EdgeBase class)
unsigned int
FASTEREdge::count_dynamic_memory() const
{
	unsigned int dynamic_memory = 0;
	dynamic_memory += two_body_energies_.size() * sizeof( core::PackerEnergy );
	dynamic_memory += parent::count_dynamic_memory();
	return dynamic_memory;
}

core::PackerEnergy FASTEREdge::get_two_body_energies_for_relaxed_states()
{
	return get_two_body_energy(
		get_faster_node( 0 )->get_relaxed_state(),
		get_faster_node( 1 )->get_relaxed_state() );
}

core::PackerEnergy
FASTEREdge::get_curr_state_energy_following_partial_state_assignment()
{
	if ( ! partial_state_assignment_ ) return curr_state_energy_;

	partial_state_assignment_ = false;
	curr_state_energy_ = get_two_body_energy(
		get_faster_node( 0 )->get_current_state(),
		get_faster_node( 1 )->get_current_state() );
	return curr_state_energy_;
}


//returns energy of perturbed state with neighbors current state
core::PackerEnergy FASTEREdge::get_energy_for_perturbed_state(
	int node,
	int nodes_perturbed_state
)
{
	return ( node == get_node_index(0) ?
		get_two_body_energy( nodes_perturbed_state, get_faster_node(1)->get_current_state() ) :
		get_two_body_energy( get_faster_node(0)->get_current_state(), nodes_perturbed_state ) );
}

void FASTEREdge::acknowledge_perturbed_state(
	int node,
	int neighbors_context_state
)
{
	int other_node = node == get_node_index(0) ? 1 : 0;
	get_faster_node( other_node )->acknowledge_neighbors_perturbed_state(
		get_edges_position_in_nodes_edge_vector( other_node ),
		neighbors_context_state);
}


void
FASTEREdge::acknowledge_participation_in_perturbation()
{
	have_contributed_deltaE_following_perturbation_ = false;
}

core::PackerEnergy
FASTEREdge::get_deltaE_for_perturbation()
{
	if ( have_contributed_deltaE_following_perturbation_ ) return core::PackerEnergy( 0.0 );
	have_contributed_deltaE_following_perturbation_ = true;
	return get_two_body_energies_for_relaxed_states() - curr_state_energy_;
}

core::PackerEnergy
FASTEREdge::get_deltaE_for_neighbor_following_perturbation( int node_index )
{
	int other_node = node_index == get_node_index(0) ? 1 : 0;
	return get_faster_node( other_node )->get_deltaE_for_relaxed_state_following_perturbation();
}


/// @details DANGER: If for some reason one were to reweight edges during simulated annealing
/// then some of the cached energies in the adjacent nodes would be out-of-date; data integrity
/// would be violated an all hell would break loose.  The same thing is true if one were to
/// change the energies on any edge during simulated annealing.  One simple solution: call
/// blanket_assign_state0 to wipe all cahced energies stored on nodes and then assign_network_state
/// to the state just before the reweighting.
/// Of course, since the annealer itself is tracking the "best" network state, you'd have to worry
/// about its data integrity as well.  General advice: don't change energies during simA.
void
FASTEREdge::set_edge_weight( Real weight )
{
	if ( weight == 0.0 ) {
		utility_exit_with_message( "Error: set edge weight to 0 not a legal operation.  Delete this edge instead" );
	}
	Real rescale = weight / edge_weight();
	two_body_energies_ *=  rescale;
	curr_state_energy_ *= rescale;
	edge_weight( weight ); // set base-class data

}


void
FASTEREdge::swap_edge_energies(
	ObjexxFCL::FArray2D< core::PackerEnergy > & new_edge_table
)
{
	if ( two_body_energies_.size1() != new_edge_table.size1() ) {
		utility_exit_with_message( "swap_edge_energies failed as size1 does not match: two_body_energies_.size1()= "
			+ utility::to_string( two_body_energies_.size1() ) + " new_edge_table.size1()= "
			+ utility::to_string( new_edge_table.size1() ) );
	}
	if ( two_body_energies_.size2() != new_edge_table.size2() ) {
		utility_exit_with_message( "swap_edge_energies failed as size2 does not match: two_body_energies_.size2()= "
			+ utility::to_string( two_body_energies_.size2() ) + " new_edge_table.size2()= "
			+ utility::to_string( new_edge_table.size2() ) );
	}
	two_body_energies_.swap( new_edge_table );
}


/// @brief main constructor: no default nor copy constructors provided.
///
/// @param num_nodes - [in] - the number of nodes in this graph
FASTERInteractionGraph::FASTERInteractionGraph(int num_nodes) :
	PrecomputedPairEnergiesInteractionGraph( num_nodes ),
	num_commits_since_last_update_(0),
	total_energy_current_state_assignment_(0),
	total_energy_alternate_state_assignment_(0),
	node_considering_alt_state_( -1 ),
	sBR_( false ),
	dBR_( false ),
	relaxed1_( -1 ),
	relaxed2_( -1 )
{}

/// @brief The FASTERIG only needs to know how many states each node has.
/// This function causes the downstream instantiation of the FASTERNodes.
void
FASTERInteractionGraph::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base )
{
	auto const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );
	for ( uint ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		set_num_states_for_node( ii, rot_sets.rotamer_set_for_moltenresidue( ii )->num_rotamers() );
	}

}

/// @brief returns the one body energy for a particular state on a node
core::PackerEnergy
FASTERInteractionGraph::get_one_body_energy_for_node_state( int node, int state)
{
	return get_faster_node( node )->get_one_body_energy( state );
}

/// @brief assigns the state of all nodes in the interaction graph to their unassigned
/// or zero states.
void FASTERInteractionGraph::blanket_assign_state_0()
{
	//a state assignment of 0 means "unassigned".
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node(ii)->assign_zero_state();
	}
	total_energy_current_state_assignment_ = 0;
	return;
}

/// @brief sets the state on node node_ind to new_state
///
/// @param node_ind - [in] - the index of the node in question
/// @param new_state - [in] - the new state the node is being assigned to
core::PackerEnergy FASTERInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	get_faster_node( node_ind )->assign_state(new_state);
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

int FASTERInteractionGraph::get_current_state_for_node( int node_ind ) const
{
	return get_faster_node( node_ind )->get_current_state();
}


/// @brief takes in a vector of states, one state per node, and sets the state for
/// each of the nodes to the specified state.
///
/// also calls "update internal energy totals" to undo any numerical noise
/// accumulated during the transition.
///
/// @param node_states - [in] - array of states, one for each node.
core::PackerEnergy FASTERInteractionGraph::set_network_state( FArray1_int & node_states)
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node( ii )->assign_state( node_states(ii) );
	}
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

/// @brief considers altering the state of a particular node; returns the
/// change in energy that the state substitution would produce
///
/// to avoid too much numerical drift from accumulating, the bookkeeping
/// arrays are updated once every 2^10 state commits
///
/// @param node_ind - [in] - the index of the node considering a state change
/// @param new_state - [in] - the new state that node is considering
/// @param alt_total_energy - [out] - the total network energy produced under the
/// new state
/// @param delta_energy - [out] - the change in energy produced under the substitution
/// @param prev_energy_for_node - [out] - the sum of the one and two body energies
///  for this node under the current state assignment
void
FASTERInteractionGraph::consider_substitution
(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
)
{
	node_considering_alt_state_ = node_ind;
	delta_energy = get_faster_node( node_ind )->
		project_deltaE_for_substitution( new_state, prev_energy_for_node );

	//numerical drift accumulates in the following assignment
	total_energy_alternate_state_assignment_ =
		total_energy_current_state_assignment_ + delta_energy;

	return;
}


/// @details to avoid too much numerical drift from accumulating, the bookkeeping
/// arrays are updated once every 2^10 state commits
core::PackerEnergy
FASTERInteractionGraph::commit_considered_substitution()
{
	get_faster_node( node_considering_alt_state_ )->
		commit_considered_substitution();
	total_energy_current_state_assignment_ =
		total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals();
	}

	return total_energy_alternate_state_assignment_;
}

core::PackerEnergy FASTERInteractionGraph::get_energy_current_state_assignment()
{
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}


/// @details Iterates across nodes and then edges to look-up the energies
/// for the current state assignmnet removing any numerical drift which
/// accumulated in the member variable total_energy_current_state_assignment_.
void FASTERInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	//std::cout << "updating internal energy totals: " << std::endl;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		//std::cout << " ig_node " << ii << " = " << ((FASTERNode *) ig_nodes_[ ii ])
		// ->get_one_body_energy_current_state();

		total_energy_current_state_assignment_ += get_faster_node( ii )->
			get_one_body_energy_current_state();
	}

	//int counter = 0;
	for ( auto iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		//std::cout << " ig_edge " << ++counter  << " =" <<
		//((FASTEREdge*) *iter)->get_current_two_body_energy();
		total_energy_current_state_assignment_ +=
			((FASTEREdge*) *iter)->get_current_two_body_energy();
	}

	//std::cout << std::endl;

	num_commits_since_last_update_ = 0;
	return;
}


int FASTERInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for ( auto iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		sum += ((FASTEREdge*) *iter)->get_two_body_table_size();
	}
	return sum;
}

/// @brief returns the amount of static memory allocated for a single FASTERInteractionGraph.
/// Does not account for any of the edges or nodes that the graph contains: that part of the
/// algorithm is implemented by the InteractionGraphBase class.
unsigned int
FASTERInteractionGraph::count_static_memory() const
{
	return sizeof( FASTERInteractionGraph );
}

/// @brief returns the amount of dynamic memory allocated for a single FASTERInteractionGraph
/// and recurses on the parent class.
/// Does not account for any of the edges or nodes that the graph contains: that part of the
/// algorithm is implemented by the InteractionGraphBase class.
unsigned int
FASTERInteractionGraph::count_dynamic_memory() const
{
	return InteractionGraphBase::count_dynamic_memory();
}

void FASTERInteractionGraph::prepare_for_FASTER()
{
	for ( auto iter = get_edge_list_begin();
			iter != get_edge_list_end();
			//note: no increment statement here
			) {
		auto next_iter = iter;
		++next_iter;
		//edges sometimes delete themselves, invalidating iterators, so
		//get the next iterator before calling prepare_for_simulated_annealing
		((FASTEREdge* ) (*iter))->prepare_for_FASTER();
		iter = next_iter;
	}

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node(ii)->prepare_for_FASTER();
	}
	return;
}

void FASTERInteractionGraph::assign_BMEC()
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node( ii )->partial_assign_state_with_lowest_one_body_energy();
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node( ii )->complete_partial_state_assignment();
	}
	update_internal_energy_totals();
}

void FASTERInteractionGraph::relax_in_current_context()
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node( ii )->relax();
	}
}

core::PackerEnergy FASTERInteractionGraph::get_energy_following_relaxation()
{
	core::PackerEnergy energy_total = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		energy_total += get_faster_node( ii )->get_one_body_energy_for_relaxed_state();
	}

	for ( auto edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		energy_total += ((FASTEREdge*) (*edge_iter))->get_two_body_energies_for_relaxed_states();
	}

	return energy_total;
}

void FASTERInteractionGraph::reject_perturbation()
{
	if ( sBR_ ) {
		get_faster_node( relaxed1_ )->reset_relaxed_for_neighbors();
		sBR_ = false;
	}

	if ( dBR_ ) {
		get_faster_node( relaxed1_ )->reset_relaxed_for_neighbors();
		get_faster_node( relaxed2_ )->reset_relaxed_for_neighbors();
		dBR_ = false;
	}
}


void FASTERInteractionGraph::commit_relaxation()
{
	sBR_ = dBR_= false;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_faster_node( ii )->partial_assign_relaxed_state(1.0);
	}

	for ( int ii = 1; ii <= get_num_nodes() ; ++ii ) {
		get_faster_node( ii )->complete_partial_state_assignment();
	}
	update_internal_energy_totals();
}

void FASTERInteractionGraph::probabilistically_commit_relaxation( Real probability )
{
	for ( int ii = 1; ii <= get_num_nodes() ; ++ii ) {
		get_faster_node( ii )->partial_assign_relaxed_state( probability );
	}
	for ( int ii = 1; ii <= get_num_nodes() ; ++ii ) {
		get_faster_node( ii )->complete_partial_state_assignment();
	}
	update_internal_energy_totals();
}

void FASTERInteractionGraph::get_current_network_state( FArray1_int & netstate )
{
	for ( int ii = 1; ii <= get_num_nodes() ; ++ii ) {
		netstate( ii ) = get_faster_node( ii )->get_current_state();
	}
}

core::PackerEnergy
FASTERInteractionGraph::perturb_sBR_and_relax(
	int node,
	int perturbed_state
)
{
	sBR_ = true;
	relaxed1_ = node;

	get_faster_node( node )->prepare_for_perturbation();
	get_faster_node( node )->set_perturbed_state( perturbed_state );
	get_faster_node( node )->tell_neighbors_to_prep_for_relaxation();
	get_faster_node( node )->relax_neighbors();

	core::PackerEnergy deltaE = get_faster_node( node )->get_deltaE_for_relaxed_state_following_perturbation();
	get_faster_node( node )->set_no_longer_perturbed();

	return deltaE;
}

core::PackerEnergy
FASTERInteractionGraph::perturb_dBR_and_relax(
	int node1,
	int perturbed_state1,
	int node2,
	int perturbed_state2
)
{
	dBR_ = true;
	relaxed1_ = node1;
	relaxed2_ = node2;

	get_faster_node( node1 )->prepare_for_perturbation();
	get_faster_node( node2 )->prepare_for_perturbation();

	get_faster_node( node1 )->set_perturbed_state( perturbed_state1 );
	get_faster_node( node2 )->set_perturbed_state( perturbed_state2 );

	get_faster_node( node1 )->tell_neighbors_to_prep_for_relaxation();
	get_faster_node( node2 )->tell_neighbors_to_prep_for_relaxation();

	get_faster_node( node1 )->relax_neighbors();
	get_faster_node( node2 )->relax_neighbors();

	core::PackerEnergy deltaE = 0;
	deltaE += get_faster_node( node1 )->get_deltaE_for_relaxed_state_following_perturbation();
	deltaE += get_faster_node( node2 )->get_deltaE_for_relaxed_state_following_perturbation();

	get_faster_node( node1 )->set_no_longer_perturbed();
	get_faster_node( node2 )->set_no_longer_perturbed();

	return deltaE;
}

int FASTERInteractionGraph::get_random_neighbor_for_node( int node )
{
	return get_faster_node( node )->get_random_neighbor();
}


void FASTERInteractionGraph::print_current_state_assignment() const
{
	std::cout << "Curr States: ";
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		std::cout << "(" << ii << ", ";
		std::cout << get_faster_node(ii)->get_current_state() << ") ";
		get_faster_node(ii)->print_internal_energies();
	}
	std::cout << std::endl;
}


/// @details For instance in a graph with 6 vertices,
/// {a,b,c,d,e,f}
/// a user may be interested in the sum of the one- and two-body energies
/// for vertices {a,b,c}.  The graph will return sum of the one body energies
/// for vertices a b and c and also any two-body energies for the edges in the
/// subgraph induced by a,b, and c.  (In this case, edges {a,b}, {a,c} and {b,c}
/// if these edges are part of the graph.  The edge {a,d} will not be counted
/// if it is part of the graph.)
/// ask the graph for the energies of the induced subgraph defined
/// by a particular group.
///
/// @param group_id - [in] - the groups for which you're interested in retrieving
///  energies of the induced subgraph
core::PackerEnergy
FASTERInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	auto esum = core::PackerEnergy( 0.0 );
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_faster_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( auto edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id ) ) {
			esum += ((FASTEREdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}
	return esum;
}


void
FASTERInteractionGraph::swap_edge_energies(
	int node1,
	int node2,
	ObjexxFCL::FArray2D< core::PackerEnergy > & new_edge_table
)
{
	get_faster_edge( node1, node2 )->swap_edge_energies( new_edge_table );
}


/// @brief factory method that instantiates a FASTERNode.
///
/// @param node_index - [in] - the index of the node being created
/// @param num_states - [in] - the total number of states for the new node
NodeBase* FASTERInteractionGraph::create_new_node( int node_index, int num_states)
{
	auto* new_node = new FASTERNode(this, node_index, num_states);
	debug_assert( new_node != nullptr );
	return new_node;
}

/// @brief factory method that instantiates a FASTEREdge
///
/// @param index1 - [in] - the smaller-indexed node this edge is incident upon
/// @param index2 - [in] - the larger-indexed node this edge is incident upon
EdgeBase* FASTERInteractionGraph::create_new_edge( int index1, int index2)
{
	return new FASTEREdge(this, index1, index2);
}

} // namespace interaction_graph
} // namespace pack
} // namespace core

