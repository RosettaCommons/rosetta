// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/DensePDInteractionGraph.cc
/// @brief  Dense Edge-Matrix, Pairwise Decomposable interaction graph class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>

// STL Headers
#include <list>
// AUTO-REMOVED #include <vector>
#include <algorithm>
#include <iostream>
// AUTO-REMOVED #include <fstream>
#include <cassert>

//Utility Headers
#include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/basic_sys_util.hh>

#include <utility/vector1.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {

//----------------------------------------------------------------------------//
//--------- Dense Pairwise Decomposable Interaction Graph Node Class ---------//
//----------------------------------------------------------------------------//

/// @brief main constructor, no default or copy constructors
///
/// allocates one-body energy array and initializes it to zero.
///
DensePDNode::DensePDNode(InteractionGraphBase * owner, int node_id, int num_states) :
	PrecomputedPairEnergiesNode( owner, node_id, num_states ),
	one_body_energies_(num_states + 1, 0.0f),
	current_state_( 0 ),
	curr_state_one_body_energy_( core::PackerEnergy( 0.0 )),
	curr_state_total_energy_( core::PackerEnergy( 0.0 )),
	alternate_state_( 0 ),
	alternate_state_one_body_energy_( core::PackerEnergy( 0.0 )),
	alternate_state_total_energy_( core::PackerEnergy( 0.0 )),
	alternate_state_is_being_considered_( false )
{}

/// @brief destructor
///
/// not responsible for any dynamically allocated memory, so node does nothing
/// it's member variables, of course, are implicitly destructed
DensePDNode::~DensePDNode()
{}

/// @brief prints a description of the node and all of it's one-body energies
void DensePDNode::print() const
{
	std::cerr << "NODE: " << get_node_index() << " with " <<
		get_num_states() << " states" << std::endl;
	for (int ii = 1; ii <= get_num_states(); ++ii)
	{
		std::cerr << "(" << ii << ", ";
		std::cerr << one_body_energies_[ ii ] << ") ";
		if ( ii % 3 == 0 ) std::cerr << std::endl;
	}
	std::cerr << std::endl  << "-----------------" << std::endl;
}


/// @brief update energy to the one-body energy for state
///
///
/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be set.
void DensePDNode::update_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
	return;
}

/// @brief set all the one-body energies for this node
///
/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void DensePDNode::update_one_body_energies( FArray1< core::PackerEnergy > & energies )
{
	assert( energies.size() == (unsigned int) get_num_states() );
	for (int ii = 1; ii <= get_num_states(); ++ii)
	{
		one_body_energies_[ ii ] = energies( ii );
	}
	return;
}

/// @brief adds energy to the one-body energy for state state
///
/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be added.
void DensePDNode::add_to_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] += energy;
	return;
}

/// @brief adds all the energies in energies to the one-body energies for this node
///
/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void DensePDNode::add_to_one_body_energies( FArray1< core::PackerEnergy > & energies )
{
	assert( energies.size() == (unsigned int) get_num_states() );
	for (int ii = 1; ii <= get_num_states(); ++ii)
	{
		one_body_energies_[ ii ] += energies( ii );
	}
	return;
}

/// @brief sets all of the one-body energies for this node to zero
void DensePDNode::zero_one_body_energies()
{
	for (int ii = 1; ii <= get_num_states(); ++ii)
	{
		one_body_energies_[ ii ] = 0;
	}
}

/// @brief returns the one body energy for a state
///
/// @param state - [in]
core::PackerEnergy DensePDNode::get_one_body_energy( int state )
{
	return one_body_energies_[ state ];
}

/// @brief prepares node for simulated annealing
///
/// updates internal edge vector + other vectorized edge information
void DensePDNode::prepare_for_simulated_annealing()
{
	if (! get_edge_vector_up_to_date() ) update_internal_vectors();
	return;
}

/// @brief assigns node's state to it's zero, or "unassigned" state.
///
/// zeros the edge-energy array, informs neighbors that it's in its unassigned
/// state
void DensePDNode::assign_zero_state()
{

	//std::cerr << "assign_state: node -  " << get_node_index() <<
	//	" new state " << 0 << "...";

	current_state_ = 0;
	alternate_state_ = 0;
	alternate_state_is_being_considered_ = false;

	curr_state_one_body_energy_ = 0.0f;
	//fills from [1] to end
	std::vector< core::PackerEnergy >::iterator position1 = curr_state_two_body_energies_.begin();
	++position1;
	std::fill( position1,
		curr_state_two_body_energies_.end(),
		0.0f);
	curr_state_total_energy_ = 0.0f;

	for (int ii = 1; ii <= get_num_incident_edges(); ++ii )
	{
		get_incident_dpd_edge(ii)->
			acknowledge_state_zeroed( get_node_index() );
	}

	return;
}


/// @brief assigns node a new_state
///
/// node updates its curr_state one and two body energies
///
/// @param new_state - [in] - the new state the node should be assigned
void DensePDNode::assign_state(int new_state)
{
	assert( new_state >= 0 && new_state <= get_num_states());

	if (new_state == 0) assign_zero_state();
	else
	{
		//std::cerr << "assign_state: node -  " << get_node_index() <<
		// " new state " << new_state << "...";
		current_state_ = new_state;
		curr_state_one_body_energy_ = one_body_energies_[ current_state_ ];
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for (int ii = 1; ii <= get_num_incident_edges(); ++ii )
		{
			get_incident_dpd_edge(ii)->acknowledge_state_change(
				get_node_index(),
				current_state_,
				curr_state_two_body_energies_[ii]);

			curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
		}
		//std::cerr<< "..done" << std::endl;
	}
	return;
}

/// @brief returns the state the node is currently assigned
int DensePDNode::get_current_state() const
{
	return current_state_;
}

/// @brief returns the one body energy for the state the node is currently assigned
core::PackerEnergy DensePDNode::get_one_body_energy_current_state()
{	return curr_state_one_body_energy_;}


/// @brief tells the node that it should change its state to the last state it was
/// asked to consider (from a call to project_deltaE_for_substitution)
///
/// updates edge energy vector, iterates across neighbors having them update
/// their edge energies.  Bookkeeping recaptures performance lost by
/// leaving energy2b structure
void DensePDNode::commit_considered_substitution()
{
	assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	curr_state_one_body_energy_ = alternate_state_one_body_energy_;
	curr_state_total_energy_ = alternate_state_total_energy_;

	//copies from [1] to end
	std::vector< core::PackerEnergy >::iterator alt_position1 = alternate_state_two_body_energies_.begin();
	++alt_position1;
	std::vector< core::PackerEnergy >::iterator curr_position1 = curr_state_two_body_energies_.begin();
	++curr_position1;

	std::copy( alt_position1,
		alternate_state_two_body_energies_.end(),
		 curr_position1 );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii )
	{
		get_incident_dpd_edge(ii)->acknowledge_substitution(
			get_node_index(),
			alternate_state_two_body_energies_[ii],
			current_state_
		);
	}

	alternate_state_is_being_considered_ = false;
	return;
}


/// @brief updates bookkeeping arrays that correspond to edge-list.
///
/// calls base class update_edge_vector function, and then proceeds to create
/// appropriate bookkeeping arrays used in simulated annealing
void DensePDNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	//aa_offsets_for_this_lookup_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_.resize( get_num_incident_edges() + 1);

	edge_matrix_ptrs_.clear();
	edge_matrix_ptrs_.reserve( get_num_incident_edges() + 1);
	edge_matrix_ptrs_.push_back( FArray2A< core::PackerEnergy >() ); //occupy the 0th position

	for (int ii = 1; ii <= get_num_incident_edges(); ++ii)
	{
		edge_matrix_ptrs_.push_back( get_incident_dpd_edge(ii)->get_edge_table_ptr() );
	}

	curr_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	alternate_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	return;
}

/// @brief outputs to standard error the bookkeeping energies for the node in its
/// current state assignment
void DensePDNode::print_internal_energies() const
{
	std::cerr << "curr_state " << current_state_ << " ";
	std::cerr << "curr_state_one_body_energy_ ";
	std::cerr << curr_state_one_body_energy_ << " ";
	std::cerr << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii)
	{
		std::cerr << "(" << get_index_of_adjacent_node( ii ) << ": " <<
			curr_state_two_body_energies_[ ii ] << ") ";
	}
	std::cerr << std::endl;
}

/// @brief removes numerical drift long stretches of efficient bookkeeping
/// produces
void DensePDNode::update_internal_energy_sums()
{
	assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = 0;
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii)
	{
		curr_state_total_energy_ +=
			get_incident_dpd_edge(ii)->get_current_two_body_energy();
	}
	curr_state_total_energy_ += curr_state_one_body_energy_;
	return;
}

/// @brief If DensePDNode is the most-derived class being used, then this function
/// will be called and will return the amount of memory statically allocated by
/// a single DensePDNode.
unsigned int
DensePDNode::count_static_memory() const
{
	return sizeof( DensePDNode );
}

/// @brief Called either by the IGBase if the DensePDNode is the most-derived class, or
/// called recursively by a derived class.  Called to account for the dynamically allocated
/// memory that this node uses.
unsigned int
DensePDNode::count_dynamic_memory() const
{
	unsigned int dynamic_memory = 0;
	dynamic_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );
	dynamic_memory += neighbors_curr_state_.size() * sizeof( int );
	dynamic_memory += edge_matrix_ptrs_.size() * sizeof ( FArray2A< core::PackerEnergy > );
	dynamic_memory += curr_state_two_body_energies_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += alternate_state_two_body_energies_.size() * sizeof ( core::PackerEnergy );
	dynamic_memory += NodeBase::count_dynamic_memory();
	return dynamic_memory;
}

//-------------------------------------------------------------------------------//
//--------- Dense Pairwise Decomposable Interaction Graph Edge Class ------------//
//-------------------------------------------------------------------------------//

/// @brief main constructor - no default nor copy constructors provided
///
/// @param owner - [in] - pointer to the graph that created this node
/// @param first_node_ind - [in] - the index of the smaller-indexed node
/// @param second_node_ind - [in] - the index of the larger-indexed node
DensePDEdge::DensePDEdge
(	InteractionGraphBase* owner,
	int first_node_ind,
	int second_node_ind
) :
	PrecomputedPairEnergiesEdge( owner, first_node_ind, second_node_ind),
	two_body_energies_(
		get_dpd_node(1)->get_num_states(),
		get_dpd_node(0)->get_num_states(),
		0.0f
	),
	curr_state_energy_( core::PackerEnergy( 0.0 ) ),
	energies_updated_since_last_prep_for_simA_( true )
{
}

/// @brief destructor.  All dynamically allocated memory is managed by the objects contained
/// inside the DensePDEdge, so there is no work to be (explicitly) done.
DensePDEdge::~DensePDEdge()
{}

/// @brief adds the input energy to the two body energy for state1 on the node with the
/// smaller index and state2 on the node with the larger index.
void DensePDEdge::add_to_two_body_energy
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
void DensePDEdge::add_to_two_body_energies
(
	FArray2< core::PackerEnergy > const & res_res_energy_array
)
{
	assert( res_res_energy_array.size1() == two_body_energies_.size1() );
	assert( res_res_energy_array.size2() == two_body_energies_.size2() );
	for ( Size ii = 1, iie = two_body_energies_.size1(); ii <= iie; ++ii ) {
		for ( Size jj = 1, jje = two_body_energies_.size2(); jj <= jje; ++jj ) {
			two_body_energies_( ii, jj ) += edge_weight() * res_res_energy_array( ii, jj );
		}
	}
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

bool DensePDEdge::get_sparse_aa_info( int, int) const {return true;} //"all amino acids are neighbors"


/// @brief Sets the two-body energy for a pair of states.  That is, it overwrites
/// whatever two-body energy there was previously for that state pair with
/// a new energy.  Ignores non-neighboring state pairs.
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
/// @param energy - [in] - the energy which replaces the old two-body energy
void DensePDEdge::set_two_body_energy
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
void DensePDEdge::clear_two_body_energy
(
	int const state1,
	int const state2
)
{
	two_body_energies_(state2,state1) = 0.0f;
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief returns the two body energy for a pair of states: 0 if those states are
/// not neighbors
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
core::PackerEnergy DensePDEdge::get_two_body_energy( int const state1, int const state2) const
{
	return two_body_energies_(state2, state1);
}

/// @brief If all of the energies for an edge have been added in, then declare the edge energies
/// final.  This may mean that the edge deletes itself.
void DensePDEdge::declare_energies_final()
{
	prepare_for_simulated_annealing();
}

/// @brief looks at all pair energies, and if they are all 0, deletes itself
void DensePDEdge::prepare_for_simulated_annealing()
{
	if ( ! energies_updated_since_last_prep_for_simA_ ) return;


	energies_updated_since_last_prep_for_simA_ = false;

	bool any_non_zero = false;
	unsigned int const num_energies = two_body_energies_.size();
	for (unsigned int ii = 0; ii < num_energies; ++ii)
	{
		if ( two_body_energies_[ ii ] != 0.0f ) { any_non_zero = true; break;}
	}

	if ( ! any_non_zero ) delete this;
}

/// @brief returns the two body energy corresponding to the current states assigned to
/// the nodes this edge is incident upon.
core::PackerEnergy DensePDEdge::get_current_two_body_energy()
{
	return curr_state_energy_;
}

/// @brief updates bookkeeping information when one of the two nodes changes its state
///
/// @param node_ind - [in] - the index of the node that changed its state
/// @param node_state - [in] - the index of the new state it assumed
/// @param new_energy - [out] - the two body energy produced  by the new state and
/// 	the current state on the other node
void
DensePDEdge::acknowledge_state_change
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
		get_dpd_node( node_not_substituted )->get_current_state();

	bool one_node_in_zero_state =
		( nodes_curr_states[0] == 0 || nodes_curr_states[1] == 0 );

	if (  one_node_in_zero_state ) {
		curr_state_energy_ = 0;
	} else {
		curr_state_energy_ = two_body_energies_( nodes_curr_states[ 1 ], nodes_curr_states[ 0 ] );
	}
	new_energy = curr_state_energy_;

	get_dpd_node( node_not_substituted )->acknowledge_neighbors_state_substitution (
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
void DensePDEdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted = ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;

	get_dpd_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		0
	);
}




/// @brief Returns a reference to the first element in the dense two-body energy
/// table.  Used to create a proxy array on the nodes for cache efficiency.
FArray2A< core::PackerEnergy > DensePDEdge::get_edge_table_ptr() {
	return FArray2A< core::PackerEnergy >( two_body_energies_( 1, 1 ),
		get_num_states_for_node( 1 ),
		get_num_states_for_node( 0 ) );
}


/// @brief returns the memory usage of the two body energy table for this edge
int DensePDEdge::get_two_body_table_size() const
{
	return two_body_energies_.size();
}

/// @brief returns sizeof DensePDEdge if this is the most-derived instance of the class
unsigned int
DensePDEdge::count_static_memory() const
{
	return sizeof( DensePDEdge );
}

/// @brief returns the amount of memory dynamically allocated by the edge and
/// recurses on its parent (in this case, the EdgeBase class)
unsigned int
DensePDEdge::count_dynamic_memory() const
{
	unsigned int dynamic_memory = 0;
	dynamic_memory += two_body_energies_.size() * sizeof( core::PackerEnergy );
	dynamic_memory += EdgeBase::count_dynamic_memory();
	return dynamic_memory;
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
DensePDEdge::set_edge_weight( Real weight )
{
	if ( weight == 0.0 ) {
		utility_exit_with_message( "Error: set edge weight to 0 not a legal operation.  Delete this edge instead" );
	}
	Real rescale = weight / edge_weight();
	two_body_energies_ *=  rescale;
	curr_state_energy_ *= rescale;
	edge_weight( weight ); // set base-class data

}

ObjexxFCL::FArray2D< core::PackerEnergy >
DensePDEdge::get_aa_submatrix_energies() const
{
	return two_body_energies_;
}

void
DensePDEdge::swap_edge_energies(
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


//----------------------------------------------------------------------------//
//--------- Dense Pairwise Decomposable Interaction Graph Class --------------//
//----------------------------------------------------------------------------//


/// @brief main constructor: no default nor copy constructors provided.
///
/// @param num_nodes - [in] - the number of nodes in this graph
DensePDInteractionGraph::DensePDInteractionGraph(int num_nodes) :
	PrecomputedPairEnergiesInteractionGraph( num_nodes ),
	num_commits_since_last_update_(0),
	total_energy_current_state_assignment_(0),
	total_energy_alternate_state_assignment_(0),
	node_considering_alt_state_( -1 )
{}

/// @brief The DensePDIG only needs to know how many states each node has.
/// This function causes the downstream instantiation of the DensePDNodes.
void
DensePDInteractionGraph::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base )
{
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );
	for ( uint ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		set_num_states_for_node( ii, rot_sets.rotamer_set_for_moltenresidue( ii )->num_rotamers() );
	}

}

/// @brief returns the one body energy for a particular state on a node
core::PackerEnergy
DensePDInteractionGraph::get_one_body_energy_for_node_state( int node, int state)
{
	return get_dpd_node( node )->get_one_body_energy( state );
}

/// @brief assigns the state of all nodes in the interaction graph to their unassigned
/// or zero states.
void DensePDInteractionGraph::blanket_assign_state_0()
{
	//a state assignment of 0 means "unassigned".
	for (int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_dpd_node(ii)->assign_zero_state();
	}
	total_energy_current_state_assignment_ = 0;
	return;
}

/// @brief sets the state on node node_ind to new_state
///
/// @param node_ind - [in] - the index of the node in question
/// @param new_state - [in] - the new state the node is being assigned to
core::PackerEnergy DensePDInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	get_dpd_node( node_ind )->assign_state(new_state);
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

/// @brief takes in a vector of states, one state per node, and sets the state for
/// each of the nodes to the specified state.
///
/// also calls "update internal energy totals" to undo any numerical noise
/// accumulated during the transition.
///
/// @param node_states - [in] - array of states, one for each node.
core::PackerEnergy DensePDInteractionGraph::set_network_state( FArray1_int & node_states)
{
	for (int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_dpd_node( ii )->assign_state( node_states(ii) );
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
///	new state
/// @param delta_energy - [out] - the change in energy produced under the substitution
/// @param prev_energy_for_node - [out] - the sum of the one and two body energies
/// 	for this node under the current state assignment
void
DensePDInteractionGraph::consider_substitution
(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
)
{
	node_considering_alt_state_ = node_ind;
	delta_energy = get_dpd_node( node_ind )->
		project_deltaE_for_substitution( new_state, prev_energy_for_node );

	//numerical drift accumulates in the following assignment
	total_energy_alternate_state_assignment_ =
		total_energy_current_state_assignment_ + delta_energy;

	return;
}


/// @details to avoid too much numerical drift from accumulating, the bookkeeping
/// arrays are updated once every 2^10 state commits
core::PackerEnergy
DensePDInteractionGraph::commit_considered_substitution()
{
	get_dpd_node( node_considering_alt_state_ )->
		commit_considered_substitution();
	total_energy_current_state_assignment_ =
		total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if (num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES) {
		update_internal_energy_totals();
	}

	return total_energy_alternate_state_assignment_;
}

core::PackerEnergy DensePDInteractionGraph::get_energy_current_state_assignment()
{
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}



/// @details Iterates across nodes and then edges to look-up the energies
/// for the current state assignmnet removing any numerical drift which
/// accumulated in the member variable total_energy_current_state_assignment_.
void DensePDInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	//std::cerr << "updating internal energy totals: " << std::endl;
	for (int ii = 1; ii <= get_num_nodes(); ++ii) {
		//std::cerr << " ig_node " << ii << " = " << ((DensePDNode *) ig_nodes_[ ii ])
		//	->get_one_body_energy_current_state();

		total_energy_current_state_assignment_ += get_dpd_node( ii )->
			get_one_body_energy_current_state();
	}

	//int counter = 0;
	for (std::list<EdgeBase*>::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter) {
		//std::cerr << " ig_edge " << ++counter  << " =" <<
		//((DensePDEdge*) *iter)->get_current_two_body_energy();
		total_energy_current_state_assignment_ +=
			((DensePDEdge*) *iter)->get_current_two_body_energy();
	}

	//std::cerr << std::endl;

	num_commits_since_last_update_ = 0;
	return;
}


int DensePDInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for (std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter) {
		sum += ((DensePDEdge*) *iter)->get_two_body_table_size();
	}
	return sum;
}

/// @brief returns the amount of static memory allocated for a single DensePDInteractionGraph.
/// Does not account for any of the edges or nodes that the graph contains: that part of the
/// algorithm is implemented by the InteractionGraphBase class.
unsigned int
DensePDInteractionGraph::count_static_memory() const
{
	return sizeof( DensePDInteractionGraph );
}

/// @brief returns the amount of dynamic memory allocated for a single DensePDInteractionGraph
/// and recurses on the parent class.
/// Does not account for any of the edges or nodes that the graph contains: that part of the
/// algorithm is implemented by the InteractionGraphBase class.
unsigned int
DensePDInteractionGraph::count_dynamic_memory() const
{
	return InteractionGraphBase::count_dynamic_memory();
}

void
DensePDInteractionGraph::swap_edge_energies(
	int node1,
	int node2,
	ObjexxFCL::FArray2D< core::PackerEnergy > & new_edge_table
)
{
	get_dpd_edge( node1, node2 )->swap_edge_energies( new_edge_table );
}


void DensePDInteractionGraph::print_current_state_assignment() const
{
	std::cerr << "Curr States: ";
	for (int ii = 1; ii <= get_num_nodes(); ++ii) {
		std::cerr << "(" << ii << ", ";
		std::cerr << get_dpd_node(ii)->get_current_state() << ") ";
		get_dpd_node(ii)->print_internal_energies();
	}
	std::cerr << std::endl;
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
DensePDInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	core::PackerEnergy esum = 0;
	for (int ii = 1; ii <= get_num_nodes(); ++ii) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_dpd_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
			&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id )) {
			esum += ((DensePDEdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}
	return esum;
}

/// @details PDInteractionGraph will return aa submatrices as requested.
bool
DensePDInteractionGraph::aa_submatrix_energies_retrievable() const
{
	return true;
}

int DensePDInteractionGraph::aatype_for_node_state(
	int,
	int
) const
{
	return 1;
}

ObjexxFCL::FArray2D< core::PackerEnergy >
DensePDInteractionGraph::get_aa_submatrix_energies_for_edge(
	int node1,
	int node2,
	int,
	int
) const
{
	return get_dpd_edge( node1, node2 )->get_aa_submatrix_energies();
}



/// @brief factory method that instantiates a DensePDNode.
///
/// @param node_index - [in] - the index of the node being created
/// @param num_states - [in] - the total number of states for the new node
NodeBase* DensePDInteractionGraph::create_new_node( int node_index, int num_states)
{
	DensePDNode* new_node = new DensePDNode(this, node_index, num_states);
	assert( new_node != NULL );
	return new_node;
}

/// @brief factory method that instantiates a DensePDEdge
///
/// @param index1 - [in] - the smaller-indexed node this edge is incident upon
/// @param index2 - [in] - the larger-indexed node this edge is incident upon
EdgeBase* DensePDInteractionGraph::create_new_edge( int index1, int index2)
{
	return new DensePDEdge(this, index1, index2);
}

} // namespace interaction_graph
} // namespace pack
} // namespace core

