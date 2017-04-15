// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/PDInteractionGraph.cc
/// @brief  Pairwise Decomposable interaction graph class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/interaction_graph/PDInteractionGraph.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <utility/assert.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2D.hh>

// STL Headers
#include <list>
#include <algorithm>
#include <iostream>

// Utility Headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

static THREAD_LOCAL basic::Tracer TR("core.pack.interaction_graph.PDInteractionGraph");

//----------------------------------------------------------------------------//
//-------- Sparse Pairwise Decomposable Interaction Graph Node Class ---------//
//----------------------------------------------------------------------------//

/// @brief main constructor, no default or copy constructors
///
/// allocates one-body energy array and initializes it to zero.
/// allocates space for sparse-matrix information.
///
PDNode::PDNode(InteractionGraphBase * owner, int node_id, int num_states) :
	PrecomputedPairEnergiesNode( owner, node_id, num_states ),
	num_aa_types_( get_pdig_owner()->get_num_aatypes() ),
	num_states_for_aatype_( num_aa_types_, 0 ),
	sparse_mat_info_for_state_( num_states + 1),
	one_body_energies_(num_states + 1, 0.0f),
	current_state_( 0 ),
	curr_state_total_energy_( 0.0 ),
	alternate_state_is_being_considered_( false )
{}

/// @brief destructor
///
/// @details not responsible for any dynamically allocated memory, so node does nothing
/// it's member variables, of course, are implicitly destructed
PDNode::~PDNode()
{}

void PDNode::print() const
{
	// std::cerr << "NODE: " << get_node_index() << " with " <<
	//  get_num_states() << " states" << std::endl;
	// for ( int ii = 1; ii <= get_num_states(); ++ii ) {
	//  std::cerr << "(" << ii << ", " <<
	//   sparse_mat_info_for_state_[ii].get_aa_type() << ", ";
	//  std::cerr <<
	//   sparse_mat_info_for_state_[ii].get_state_ind_for_this_aa_type() << ", ";
	//  std::cerr << one_body_energies_[ ii ] << ") ";
	//  if ( ii % 3 == 0 ) std::cerr << std::endl;
	// }
	// std::cerr << std::endl  << "-----------------" << std::endl;

	std::cout << "Node: " << get_node_index() << " with curr_state_ " << current_state_ << " 1b: ";
	std::cout << curr_state_one_body_energy_ << std::endl;
	for ( int ii = 0; (Size) ii < curr_state_two_body_energies_.size(); ++ii ) {
		std::cout << "  2b " << get_index_of_adjacent_node( ii ) << ": " << curr_state_two_body_energies_[ ii ] << std::endl;
	}
	std::cout << " and alternate state: " << alternate_state_ << " 1b: " << alternate_state_one_body_energy_ << std::endl;
	for ( int ii = 0; (Size) ii < alternate_state_two_body_energies_.size(); ++ii ) {
		std::cout << "  2b " << get_index_of_adjacent_node( ii ) << ": " << alternate_state_two_body_energies_[ ii ] << std::endl;
	}

}

/// @details The graph doesn't know anything specific about amino acid types, for
/// instance, nothing specially distinguishes glycine and threonine.
/// The sparse-matrix representation for the PDEdge two-body energy tables
/// requires knowing when two states are of the *same* amino acid type and
/// when they are *different* amino acid types.
///
/// @param aatypes_for_state - [in] - amino acid type for each state in the node
void PDNode::set_amino_acid_types( std::vector< int > const & aatypes_for_state)
{
	debug_assert(aatypes_for_state.size() == (Size)(get_num_states() + 1) );

	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		debug_assert( aatypes_for_state[ii] > 0 && aatypes_for_state[ii] <= (int)num_states_for_aatype_.size());

		sparse_mat_info_for_state_[ii].set_aa_type( aatypes_for_state[ii] );

		++(num_states_for_aatype_[sparse_mat_info_for_state_[ii].get_aa_type() ] );

		sparse_mat_info_for_state_[ii].set_state_ind_for_this_aa_type(
			num_states_for_aatype_[sparse_mat_info_for_state_[ii].get_aa_type()]);
	}
	return;
}

int
PDNode::aatype_for_state( int state ) const
{
	return sparse_mat_info_for_state_[ state ].get_aa_type();
}


/// @details Used by AminoAcidNeighborSparseMatrix.  The FArray must persist for
/// as long as any AminoAcidNeighborSparseMatrix points to it.
///
utility::vector1< int > const &
PDNode::get_num_states_for_aa_types() const
{
	return num_states_for_aatype_;
}


/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be set.
void PDNode::update_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] = energy;
	return;
}

/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void PDNode::update_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies )
{
	debug_assert( energies.size() == (unsigned int) get_num_states() );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] = energies( ii );
	}
	return;
}


/// @param state - [in] - one-based index of the state
/// @param energy - [in] - the energy that should be added.
///
void PDNode::add_to_one_body_energy( int state, core::PackerEnergy energy )
{
	one_body_energies_[ state ] += energy;
	return;
}

/// @param energies - [in] - the array of energies. Must hold num_states_ entries
void PDNode::add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies )
{
	debug_assert( energies.size() == (unsigned int) get_num_states() );
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] += energies( ii );
	}
	return;
}

void PDNode::zero_one_body_energies()
{
	for ( int ii = 1; ii <= get_num_states(); ++ii ) {
		one_body_energies_[ ii ] = 0;
	}
}


/// @param state - [in]
core::PackerEnergy PDNode::get_one_body_energy( int state )
{
	return one_body_energies_[ state ];
}


/// @details updates internal edge vector + other vectorized edge information
void PDNode::prepare_for_simulated_annealing()
{
	if ( ! get_edge_vector_up_to_date() ) update_internal_vectors();
	return;
}

/*
/// @brief old method for figuring out how much memory the IG is using; deprecated
unsigned int
PDNode::getMemoryUsageInBytes() const
{
unsigned int total_memory = 0;
total_memory += num_states_for_aatype_.size() * sizeof( int );
total_memory += sparse_mat_info_for_state_.size() * sizeof( SparseMatrixIndex );
total_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );
total_memory += neighbors_curr_state_.size() * sizeof (int );

total_memory += aa_offsets_for_edges_.size() * sizeof( int );
total_memory += num_states_for_aa_type_for_higher_indexed_neighbor_.size() * sizeof( int );
total_memory += neighbors_curr_state_.size() * sizeof( int );
total_memory += neighbors_curr_state_sparse_info_.size() * sizeof( SparseMatrixIndex );
total_memory += edge_matrix_ptrs_.size() * sizeof( FArray1A< core::PackerEnergy > );

total_memory += curr_state_two_body_energies_.size() * sizeof( float );
total_memory += alternate_state_two_body_energies_.size() * sizeof( float );

total_memory += sizeof( PDNode );
//total_memory += NodeBase::getMemoryUsageInBytes();
return total_memory;
}
*/


/// @details zeros the edge-energy array, informs neighbors that it's in its unassigned
/// state
///
void PDNode::assign_zero_state()
{

	//std::cerr << "assign_state: node -  " << get_node_index() <<
	// " new state " << 0 << "...";

	current_state_ = 0;
	alternate_state_ = 0;
	alternate_state_is_being_considered_ = false;

	curr_state_one_body_energy_ = 0.0f;
	//fills from [1] to end
	std::vector< float >::iterator position1 = curr_state_two_body_energies_.begin();
	++position1;
	std::fill( position1,
		curr_state_two_body_energies_.end(),
		0.0f);
	curr_state_total_energy_ = 0.0f;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_pd_edge(ii)->
			acknowledge_state_zeroed( get_node_index() );
	}

	return;
}


/// @note updates its curr_state one and two body energies
///
/// @param new_state - [in] - the new state the node should be assigned
///
void PDNode::assign_state(int new_state)
{
	debug_assert( new_state >= 0 && new_state <= get_num_states());

	if ( new_state == 0 ) {
		assign_zero_state();
	} else {
		//std::cerr << "assign_state: node -  " << get_node_index() <<
		// " new state " << new_state << "...";
		current_state_ = new_state;
		curr_state_sparse_mat_info_ =
			sparse_mat_info_for_state_[ current_state_ ];
		curr_state_one_body_energy_ = one_body_energies_[ current_state_ ];
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_pd_edge(ii)->acknowledge_state_change(
				get_node_index(),
				current_state_,
				curr_state_sparse_mat_info_,
				curr_state_two_body_energies_[ii]);

			curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
		}
		//std::cerr<< "..done" << std::endl;
	}
	return;
}

int PDNode::get_current_state() const
{
	return current_state_;
}


float PDNode::get_one_body_energy_current_state() const
{
	return curr_state_one_body_energy_;
}


/// @brief tells the node that it should change its state to the last state it was
/// asked to consider (from a call to project_deltaE_for_substitution)
///
/// updates edge energy vector, iterates across neighbors having them update
/// their edge energies.  Bookkeeping recaptures performance lost by
/// leaving energy2b structure
void PDNode::commit_considered_substitution()
{
	debug_assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	curr_state_sparse_mat_info_ = alt_state_sparse_mat_info_;
	curr_state_one_body_energy_ = alternate_state_one_body_energy_;
	curr_state_total_energy_ = alternate_state_total_energy_;

	//copies from [1] to end
	std::vector< float >::iterator alt_position1 = alternate_state_two_body_energies_.begin();
	++alt_position1;
	std::vector< float >::iterator curr_position1 = curr_state_two_body_energies_.begin();
	++curr_position1;

	std::copy( alt_position1,
		alternate_state_two_body_energies_.end(),
		curr_position1 );

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_pd_edge(ii)->acknowledge_substitution(
			get_node_index(),
			alternate_state_two_body_energies_[ii],
			current_state_,
			curr_state_sparse_mat_info_
		);
	}

	alternate_state_is_being_considered_ = false;
	return;
}


/// @brief updates bookkeeping arrays that correspond to edge-list.
///
/// calls base class update_edge_vector function, and then proceeds to create
/// appropriate bookkeeping arrays used in simulated annealing
///
/// It's possible that a derived class created an Edge, but that Edge doesn't have any two body energies. For example,
/// in the HPatchIG, two residues (Nodes) may have SASA overlap, but be outside of the 5.5A energy function cutoff.
/// To handle this case, check the size of the two-body energy table before trying to dereference it. (ronj)
///
void PDNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	//aa_offsets_for_this_lookup_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_sparse_info_.resize( get_num_incident_edges() + 1);

	edge_matrix_ptrs_.clear();
	edge_matrix_ptrs_.reserve( get_num_incident_edges() + 1);
	edge_matrix_ptrs_.push_back( ObjexxFCL::FArray1A< core::PackerEnergy >() ); //occupy the 0th position

	aa_offsets_for_edges_.dimension(
		num_aa_types_, get_num_incident_edges(), num_aa_types_);
	num_states_for_aa_type_for_higher_indexed_neighbor_.dimension(
		num_aa_types_, get_num_edges_to_larger_indexed_nodes());

	//copy offsets from edges
	//int neighb_aa_offset =
	// num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1);
	int count_neighbs_with_higher_indices = 0;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		neighbors_curr_state_sparse_info_[ii].set_aa_type( 1 );

		// Edge::get_edge_table_ptr() calls getMatrixPointer() on the AminoAcidNeighborSparseMatrix instance two_body_energies_
		// kept on the Edge. getMatrixPointer() returns a reference to the first element in this table. Problem is that
		// if there are no two-body energies, dereferencing the pointer that's kept in AANSM causes an FArray operator()
		// out-of-bounds error. Instead, reorder these lines so that we first check the two-body table size, and only call
		// get_edge_table_ptr() if the two-body energy table is nonzero.  If it's all zeros, just add an FArray1A object
		// that's been default constructed (the default constructor just sets the internal pointers to NULL) to the
		// edge_matrix_ptrs_. (ronj)
		int edge_table_size = get_incident_pd_edge(ii)->get_two_body_table_size();
		if ( edge_table_size != 0 ) {
			float & edge_table_ref = get_incident_pd_edge(ii)->get_edge_table_ptr();
			edge_matrix_ptrs_.push_back( ObjexxFCL::FArray1A< core::PackerEnergy >( edge_table_ref ));
			edge_matrix_ptrs_[ii].dimension( edge_table_size );
		} else {
			edge_matrix_ptrs_.push_back( ObjexxFCL::FArray1A< core::PackerEnergy >() );
		}

		ObjexxFCL::FArray2D_int const & edge_aa_neighb_offsets =
			get_incident_pd_edge(ii)->get_offsets_for_aatypes();
		utility::vector1< int > const & neighb_num_states_per_aa =
			get_incident_pd_edge(ii)->get_second_node_num_states_per_aa();

		if ( get_node_index() < get_index_of_adjacent_node(ii) ) {
			++count_neighbs_with_higher_indices;
			for ( int jj = 1; jj <= num_aa_types_; ++jj ) {
				for ( int kk = 1; kk <= num_aa_types_; ++kk ) {
					aa_offsets_for_edges_(kk, ii, jj) = edge_aa_neighb_offsets(kk, jj);
				}
				num_states_for_aa_type_for_higher_indexed_neighbor_(
					jj, count_neighbs_with_higher_indices) =
					neighb_num_states_per_aa[ jj ];
				//++neighb_aa_offset;
			}
		} else {
			for ( int jj = 1; jj <= num_aa_types_; ++jj ) {
				for ( int kk = 1; kk <= num_aa_types_; ++kk ) {
					aa_offsets_for_edges_(kk, ii, jj) =
						edge_aa_neighb_offsets(jj, kk);
				}
			}
		}
	}

	curr_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	alternate_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	return;
}

/// @brief - allow derived class to "drive" through the deltaE calculation
void
PDNode::calc_deltaEpd( int alternate_state )
{
	//std::cout << "PDNode::calc_deltaEpd" << std::endl;

	float dummy(0.0f);
	project_deltaE_for_substitution( alternate_state, dummy );
}


/// @brief returns the sparse matrix information for a paricular state
///
/// @param state - [in] - the state whose information is requested
SparseMatrixIndex const &
PDNode::get_sparse_mat_info_for_state(int state) const
{
	debug_assert( state > 0 && state <= get_num_states());
	return sparse_mat_info_for_state_[ state ];
}

/// @brief returns the sparse matrix information for the current state
///
SparseMatrixIndex const &
PDNode::get_sparse_mat_info_for_curr_state() const
{
	return curr_state_sparse_mat_info_;
}

/// @brief returns the number of states that are of a particular amino acid type
///
/// @param aa_type - [in] - the amino acid type in question
///
int
PDNode::get_num_states_for_aa_type(int aa_type) const
{
	return num_states_for_aatype_[ aa_type ];
}

/// @brief outputs to standard error the bookkeeping energies for the node in its
/// current state assignment
///
void PDNode::print_internal_energies() const
{
	std::cerr << "curr_state " << current_state_ << " ";
	std::cerr << "curr_state_sparse_mat_info_ ";
	std::cerr << curr_state_sparse_mat_info_.get_aa_type() << " ";
	std::cerr << curr_state_sparse_mat_info_.get_state_ind_for_this_aa_type() << " ";
	std::cerr << "curr_state_one_body_energy_ ";
	std::cerr << curr_state_one_body_energy_ << " ";
	std::cerr << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		std::cerr << "(" << get_index_of_adjacent_node( ii ) << ": " <<
			curr_state_two_body_energies_[ ii ] << ") ";
	}
	std::cerr << std::endl;
}

/// @brief removes numerical drift long stretches of efficient bookkeeping
/// produces
///
void PDNode::update_internal_energy_sums()
{
	debug_assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = 0;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_total_energy_ +=
			get_incident_pd_edge(ii)->get_current_two_body_energy();
	}
	curr_state_total_energy_ += curr_state_one_body_energy_;
	return;
}

unsigned int
PDNode::count_static_memory() const
{
	return sizeof ( PDNode );
}

unsigned int
PDNode::count_dynamic_memory() const
{
	unsigned int total_memory = 0;
	total_memory += num_states_for_aatype_.size() * sizeof( int );
	total_memory += sparse_mat_info_for_state_.size() * sizeof( SparseMatrixIndex );
	total_memory += one_body_energies_.size() * sizeof( core::PackerEnergy );
	total_memory += neighbors_curr_state_.size() * sizeof (int );

	total_memory += aa_offsets_for_edges_.size() * sizeof( int );
	total_memory += num_states_for_aa_type_for_higher_indexed_neighbor_.size() * sizeof( int );
	total_memory += neighbors_curr_state_.size() * sizeof( int );
	total_memory += neighbors_curr_state_sparse_info_.size() * sizeof( SparseMatrixIndex );
	total_memory += edge_matrix_ptrs_.size() * sizeof( ObjexxFCL::FArray1A< core::PackerEnergy > );

	total_memory += curr_state_two_body_energies_.size() * sizeof( float );
	total_memory += alternate_state_two_body_energies_.size() * sizeof( float );

	total_memory += NodeBase::count_dynamic_memory();
	return total_memory;

}


/*

/// @brief ready a node for the output-to-file process
///
/// allocates space for an additional array that counts the number of states
/// for each amino acid
///
void PDNode::prepare_to_write_to_file()
{
initialize_aa_for_state_array();
}


/// @brief counts the number of states for each amion acid
void PDNode::initialize_aa_for_state_array()
{
aa_types_for_instance_states_.dimension( get_num_states() );
for (int ii = 1; ii <= get_num_states(); ++ii)
{
aa_types_for_instance_states_(ii) = sparse_mat_info_for_state_[ii].get_aa_type();
}
}

/// @brief deallocates extra memory allocated before writing to a file
void PDNode::clean_up_after_writing_to_file()
{
aa_types_for_instance_states_.dimension(0);
}

/// @brief allocates extra space required of reading to a file
///
/// this graph instance may not correspond perfectly to the rotamers
/// described in the file.  Before reading pair energies from the file
/// a correspondnece must be found between the "instance states" and
/// and the "file states".  RotamerSets.cc finds the correspondence.
/// Each node keeps track of the correspondence - which instance states
/// correspond to which file states, and which instance states have no
/// corresponding file states.  A vertex also keeps track of which
/// amino acid type for each of file state.  Amino acid type information is
/// necessary for reading in the sparse table that relies on the amino-acid
/// neighbor optimization.
///
/// @param num_states_for_node_in_file - [in] - the number of states on the corresponding
///   node from the file
///
void PDNode::prepare_to_read_energies_from_file(int num_states_for_node_in_file)
{
int num_aa_types_in_file =
get_pdig_owner()->get_num_file_aatypes();

num_states_in_file_ = num_states_for_node_in_file;
instance_states_2_file_states_.dimension( get_num_states() );
file_states_2_instance_states_.dimension( num_states_for_node_in_file );
aa_types_for_file_states_.dimension( num_states_for_node_in_file );
num_file_states_for_aa_.dimension( num_aa_types_in_file );

instance_states_2_file_states_ = -1;
file_states_2_instance_states_ = -1;
aa_types_for_file_states_ = -1;
num_file_states_for_aa_ = 0;

initialize_aa_for_state_array(); //apl borrowing functionality from wri
}

/// @brief deallocate extra tables after having written to a file
void PDNode::clean_up_after_reading_energies_from_file()
{
instance_states_2_file_states_.dimension(0);
file_states_2_instance_states_.dimension(0);
aa_types_for_file_states_.dimension(0);
aa_types_for_instance_states_.dimension(0);
num_file_states_for_aa_.dimension(0);
}

/// @brief set the amino acid type for a file-state
///
/// This method assumes that each file state will have its aa type set
/// before edge energies are read.  The structure of the interaction
/// graph file ensures that the amino acid type information is present for
/// each file state.  During RotamerSets.cc::read_rotamer_sets_from_file
/// the information must be transmitted to the interaction graph.
/// If this code quit rosetta, it's because the intereaction graph file
/// holds an amino-acid type out of bounds.
///
/// @param file_state - [in] - the index of the file state
/// @param aa - [in] - the amino acid type of the file state
///
void PDNode::set_aa_for_file_state(int file_state, int aa )
{
if ( (unsigned int) aa > num_file_states_for_aa_.size() || aa <= 0)
{
std::cerr << "Error in interaction graph file: amino acid type out";
std::cerr << " of range on node: " << get_node_index() << " for file state ";
std::cerr <<  file_state << ": aa = " << aa << std::endl;
utility_exit();
}
aa_types_for_file_states_( file_state ) = aa;
++num_file_states_for_aa_( aa );
}

/// @brief declares that an instance state corresponds to a file state.
///
/// Trouble arises if rotamers are repeated, and the code that detects the
/// correspondence doesn't prevent multiple correspondences (multiple file
/// states to a single instance state or multiple instance states to a
/// single file state.)  The way I've written the RotamerSets.cc::
/// find_correspondence_between_file_and_instance() code enforces the bijection.
///
/// @param instance_state - [in] - the index of the instance state
/// @param state_from_file - [in] - the index of the file state
///
void PDNode::set_instance_state_correspondence
(
int instance_state,
int state_from_file
)
{
//apl Enforce bijection:
//apl Much easier to handle I/O when no to file-states map to a single
//apl instance state, and no two instance states map to a single file state.
//apl If a protocol creates two identical rotamers before reading the
//apl interaction graph file, then that protocol must also
//apl have generated the same identical rotamers before writing the
//apl interaction graph file.
if ( instance_states_2_file_states_( instance_state ) != -1
|| file_states_2_instance_states_( state_from_file ) != -1 )
{
std::cerr << "Reading Interaction Graph from File: Bijection Failure" << std::endl;
std::cerr << "Node: " << get_node_index() << " instance_state " << instance_state;
std::cerr << "file_state: " << state_from_file << std::endl;
std::cerr << "First Correspondence: instance_states_2_file_states_( ";
std::cerr << instance_state << " ) = " << instance_states_2_file_states_( instance_state );
std::cerr << std::endl << "First Correspondence: file_states_2_instance_staets_( ";
std::cerr << state_from_file << " ) = " << file_states_2_instance_states_( state_from_file );
std::cerr << std::endl;
utility_exit();
}

instance_states_2_file_states_( instance_state ) = state_from_file;
file_states_2_instance_states_( state_from_file ) = instance_state;
}


/// @brief return the file state that corresponds to an instance state
///
/// @param instance_state - [in] - the index of the instance state
///
int PDNode::get_correspondence_for_state( int instance_state )
{
return instance_states_2_file_states_( instance_state );
}

/// @brief returns the number of instance states that did not correspond to any
/// file state.  The rotamers these states correspond to were absent from
/// the file, and their pair energies must be computed.
///
int PDNode::get_num_rots_absent_from_file()
{
if ( instance_states_2_file_states_.size() == 0 )
{
//no correspondence found with input file
return get_num_states();
}

int count_absent = 0;
for (int ii = 1; ii <= get_num_states(); ++ii)
{
if (instance_states_2_file_states_(ii) == -1 )
{
++count_absent;
}
}
return count_absent;
}

/// @brief writes the index of the instance states with no matching file states into
/// the input FArray which should be large enough to hold them.
/// If the node itself had no matching node in the input file, then all
/// rotamers were absent from the file.
///
void PDNode::get_absent_rots( FArray1_int & rots_absent )
{
if ( instance_states_2_file_states_.size() == 0 ) {
//no correspondence found with input file
for (int ii = 1; ii <= get_num_states(); ++ii) {
rots_absent(ii) = ii;
}
return;
}

int count_absent = 0;
for (int ii = 1; ii <= get_num_states(); ++ii) {
if (instance_states_2_file_states_(ii) == -1 ) {
++count_absent;
rots_absent( count_absent ) = ii;
}
}
}

/// @brief returns the number of file states for this node
int PDNode::get_num_states_in_file()
{
return num_states_in_file_;
}

/// @brief returns a reference to the first entry in the aa_types_for_file_states_
/// array so that an FArray1A can be constructed.
///
/// used by a PDEdge while it's reading energies from a file
int & PDNode::get_aatypes_for_file_states()
{
return aa_types_for_file_states_(1);
}

/// @brief returns a reference to the first entry in the aa_types_for_instance_states_
/// array so that an FArray1A can be constructed.
///
/// used by a PDEdge while it's reading energies from a file
int & PDNode::get_aatypes_for_states()
{
return aa_types_for_instance_states_(1);
}

/// @brief returns a reference to the first entry in the num_file_states_for_aa_
/// array so that an FArray1A can be constructed.
///
/// used by a PDEdge while it's reading energies from a file
///
int & PDNode::get_num_file_states_for_aa()
{
return num_file_states_for_aa_(1);
}

/// @brief returns a reference to the first entry in the file_states_2_instance_states_
/// array so that an FArray1A can be constructed.
///
/// used by a PDEdge while it's reading energies from a file
///
int & PDNode::get_file_states_2_instance_states_array()
{
return file_states_2_instance_states_(1);
}

/// @brief returns true if a node corresponds to one of the nodes described in the
/// input file, and false if that node does not correspond to any.
///
bool PDNode::get_node_corresponded_to_file_node()
{
return (instance_states_2_file_states_.size() != 0 );
}

*/

//-------------------------------------------------------------------------------//
//-------- Sparse Pairwise Decomposable Interaction Graph Edge Class ------------//
//-------------------------------------------------------------------------------//

/// @brief main constructor - no default nor copy constructors provided
///
/// @param owner - [in] - pointer to the graph that created this node
/// @param first_node_ind - [in] - the index of the smaller-indexed node
/// @param second_node_ind - [in] - the index of the larger-indexed node
///
PDEdge::PDEdge
( InteractionGraphBase* owner,
	int first_node_ind,
	int second_node_ind
) :
	PrecomputedPairEnergiesEdge( owner, first_node_ind, second_node_ind),
	//energy_table_size_(0)
	two_body_energies_(
	get_pd_node(0)->get_num_states_for_aa_types(),
	get_pd_node(1)->get_num_states_for_aa_types()
	),
	energies_updated_since_last_prep_for_simA_( true )
{
	force_all_aa_neighbors();
}

/// @brief destructor
PDEdge::~PDEdge()
{}

/// @brief allocates two-body energy table based on amino-acid neighbor relationships
/// and initializes the table to 0.
///
/// @param sparse_conn_info - [in] - a MAX_AA x MAX_AA 2D array where each "true" entry
///  means that the corresponding amino acid pair are neighbors.
///
/// @remarks idea borrowed from energy2b implementation
///
void PDEdge::set_sparse_aa_info(ObjexxFCL::FArray2_bool const & sparse_conn_info)
{
	two_body_energies_.set_sparse_aa_info( sparse_conn_info );
	energies_updated_since_last_prep_for_simA_ = true;
}

/// @brief returns whether two amino acid types are represented as neighbors
bool PDEdge::get_sparse_aa_info( int node1aa, int node2aa ) const
{
	return two_body_energies_.get_sparse_aa_info( node1aa, node2aa );
}

/// @brief re-allocates two-body energy table after forcing a pair of amino acids
/// to become neighbors that were not initially declared to be neighbors
///
/// @param node1aa - [in] - the amino acid type for the node with the smaller index
/// @param node2aa - [in] - the amino acid type for the node with the larger index
///
void PDEdge::force_aa_neighbors(int node1aa, int node2aa)
{
	two_body_energies_.force_aa_neighbors( node1aa, node2aa );
	energies_updated_since_last_prep_for_simA_ = true;
}

/// @brief re-allocates two-body energy table after forcing a pair of amino acids
/// to become neighbors that were not initially declared to be neighbors
///
/// @param node1aa - [in] - the amino acid type for the node with the smaller index
/// @param node2aa - [in] - the amino acid type for the node with the larger index
///
void PDEdge::force_all_aa_neighbors()
{
	two_body_energies_.force_all_aa_neighbors();
	energies_updated_since_last_prep_for_simA_ = true;
}

/// @brief adds the input energy to the two body energy for state1 on the node with the
/// smaller index and state2 on the node with the larger index so long as
/// the amion acid types of those states have been previously declared amino
/// acid neighbors.  Any energies for non-neighboring states are ignored.
///
void PDEdge::add_to_two_body_energy
(
	int const state1,
	int const state2,
	float const energy
)
{
	two_body_energies_.add(
		get_pd_node(0)->get_sparse_mat_info_for_state(state1),
		get_pd_node(1)->get_sparse_mat_info_for_state(state2),
		edge_weight() * energy);
	energies_updated_since_last_prep_for_simA_ = true;
	return;
}

/// @brief Adds all the energies stored in the oversized_res_res_energy array to the
/// two body energy table for those states whose amion acid types were
/// previoudsly declared to be amino-acid neighbors.  The res-res array
/// should have the dimension (node1->get_num_states() x node2->get_num_states());
///
/// @param res_res_energy_array - [in] - an array containing the state pair
///  energies
///
void PDEdge::add_to_two_body_energies
(
	ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array
)
{
	for ( int ii = 1; ii <= get_num_states_for_node(0); ++ii ) {
		SparseMatrixIndex const & state1_sparse_info = get_pd_node(0)
			->get_sparse_mat_info_for_state( ii );
		for ( int jj = 1; jj <= get_num_states_for_node(1); ++jj ) {
			SparseMatrixIndex const & state2_sparse_info = get_pd_node(1)
				->get_sparse_mat_info_for_state( jj );

			two_body_energies_.add(
				state1_sparse_info, state2_sparse_info,
				edge_weight() * res_res_energy_array( jj, ii )
			);
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
///
void PDEdge::set_two_body_energy
(
	int const state1,
	int const state2,
	float const energy
)
{
	two_body_energies_.set(
		get_pd_node(0)->get_sparse_mat_info_for_state(state1),
		get_pd_node(1)->get_sparse_mat_info_for_state(state2),
		edge_weight() * energy
	);
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
///
void PDEdge::clear_two_body_energy
(
	int const state1,
	int const state2
)
{
	two_body_energies_.set(
		get_pd_node(0)->get_sparse_mat_info_for_state(state1),
		get_pd_node(1)->get_sparse_mat_info_for_state(state2),
		0.
	);

	return;
}

/// @brief returns the two body energy for a pair of states: 0 if those states are
/// not neighbors
///
/// @param state1 - [in] - state index for the node with the smaller index
/// @param state2 - [in] - state index for the node with the larger index
///
float PDEdge::get_two_body_energy( int const state1, int const state2) const
{
	return two_body_energies_.get(
		get_pd_node(0)->get_sparse_mat_info_for_state(state1),
		get_pd_node(1)->get_sparse_mat_info_for_state(state2));
}

/// @brief When all the energies that are going to be stored in an edge have been
/// placed in it, the edge may save some memory by shrinking its AminoAcidNeighborSparseMatrix.
/// This method instructs the edge to do so.
void PDEdge::declare_energies_final()
{
	prepare_for_simulated_annealing();
}

/// @brief reduces the size of the pair-energy table if any amino-acid-neighbor
/// submatrices hold nothing but 0's
///
/// since the drop_zero_submatrices_where_possible() method of the AANSM is
/// somewhat time consuming, and since it can only reduce memory use / simA
/// running time on the first execution following an update to the two-body
/// energies, the PDEdge member variable energies_updated_since_last_prep_
/// for_simA ensures that the AANSM method is only called once following
/// the update of any RPEs.
///
void PDEdge::prepare_for_simulated_annealing()
{
	prepare_for_simulated_annealing_no_deletion();
	if ( two_body_energies_.get_table_size() == 0 ) delete this;
}

/*
unsigned int
PDEdge::getMemoryUsageInBytes() const
{
unsigned int total_memory = 0;
total_memory += two_body_energies_.get_table_size() * sizeof( int );
total_memory += two_body_energies_.get_offset_table_size_in_bytes();
total_memory += sizeof( PDEdge );
return total_memory;
}
*/

/// @brief returns the two body energy corresponding to the current states assigned to
/// the nodes this edge is incident upon.
///
float PDEdge::get_current_two_body_energy()
{
	return curr_state_energy_;
}

/// @brief updates bookkeeping information when one of the two nodes changes its state
///
/// @param node_ind - [in] - the index of the node that changed its state
/// @param node_state - [in] - the index of the new state it assumed
/// @param new_state_spare_info - [in] - the sparse-matrix information for the state
/// @param new_energy - [out] - the two body energy produced  by the new state and
///  the current state on the other node
///
void
PDEdge::acknowledge_state_change
(
	int node_ind,
	int new_state,
	SparseMatrixIndex const & new_state_sparse_info,
	float & new_energy
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	int nodes_curr_states[2];
	SparseMatrixIndex nodes_curr_states_sparse_info[2];

	nodes_curr_states[ node_substituted ] = new_state;
	nodes_curr_states_sparse_info[ node_substituted ] = new_state_sparse_info;

	nodes_curr_states[ node_not_substituted ] =
		get_pd_node( node_not_substituted )->get_current_state();
	nodes_curr_states_sparse_info[ node_not_substituted ] =
		get_pd_node( node_not_substituted )->
		get_sparse_mat_info_for_curr_state();

	bool one_node_in_zero_state =
		( nodes_curr_states[0] == 0 || nodes_curr_states[1] == 0 );

	if (  one_node_in_zero_state ) {
		curr_state_energy_ = 0;
	} else {
		curr_state_energy_ = two_body_energies_.get(
			nodes_curr_states_sparse_info[0],
			nodes_curr_states_sparse_info[1]);
	}
	new_energy = curr_state_energy_;

	get_pd_node( node_not_substituted )->
		acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		new_state,
		new_state_sparse_info
	);

	return;
}


/// @param node_ind - [in] - the index of the node that has just entered its 0 state
void PDEdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;
	SparseMatrixIndex dummy_sparse_info;
	dummy_sparse_info.set_aa_type( 1 );
	dummy_sparse_info.set_state_ind_for_this_aa_type(1);

	get_pd_node( node_not_substituted )->
		acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		0,
		dummy_sparse_info
	);
	return;
}


ObjexxFCL::FArray2D_int const &
PDEdge::get_offsets_for_aatypes( )
{
	return two_body_energies_.getAANeighborOffsets();
}

utility::vector1< int > const &
PDEdge::get_second_node_num_states_per_aa()
{
	return get_pd_node(1)->get_num_states_for_aa_types();
}

/// @brief Returns a reference to the first element in the sparse two-body energy
/// table.  Used to create a proxy array on the nodes for cache efficiency.
///
float & PDEdge::get_edge_table_ptr() {
	//std::cout << "PDEdge: get_edge_table_ptr(): two_body_energies_.size(): " << two_body_energies_.get_table_size() << std::endl;
	return two_body_energies_.getMatrixPointer();
}

unsigned int
PDEdge::count_static_memory() const
{
	return sizeof ( PDEdge );
}

unsigned int
PDEdge::count_dynamic_memory() const
{
	unsigned int total_memory = 0;
	total_memory += two_body_energies_.get_table_size() * sizeof( int );
	total_memory += two_body_energies_.get_offset_table_size_in_bytes();
	total_memory += EdgeBase::count_dynamic_memory();
	return total_memory;
}

ObjexxFCL::FArray2D< core::PackerEnergy >
PDEdge::get_aa_submatrix_energies(
	int node1aa,
	int node2aa
) const
{
	return two_body_energies_.get_aa_submatrix_energies( node1aa, node2aa );
}


/// @details DANGER: If for some reason one were to reweight edges during simulated annealing
/// then some of the cached energies in the adjacent nodes would be out-of-date; data integrity
/// would be violated an all hell would break loose.  The same thing is true if one were to
/// change the energies on any edge during simulated annealing.  One simple solution: call
/// blanket_assign_state0 to wipe all cahced energies stored on nodes and then assign_network_state
/// to the state just before the reweighting.
/// Of course, since the annealer itself is tracking the "best" network state, one would have to worry
/// about its data integrity as well.  General advice: don't change energies during simA.
void
PDEdge::set_edge_weight( Real weight )
{
	Real rescale = weight / edge_weight();
	two_body_energies_.scale( rescale );
	edge_weight( weight ); // set base-class data
}


/*

/// @brief Reads the energies for an edge described in the input file and
/// copies the energies that correspond to a pair of instance states into
/// the edge energy table.
///
/// Energies are represented in a binary file.  They are stored using the
/// same amino-acid neighbor technique that saves so much memory.  The
/// first part of the information in the file for this edge is the amino
/// acid neighbor connectivity information; 400 bools ( num_aa ^ 2).
/// From the number of file states for each amino acid, and the connectivity
/// the number of edge energies present can be readily computed.
/// This method allocates a large buffer to read in all the energies at once,
/// minimizing the number of disk reads.  It then walks through the buffer
/// and writes the appropriate energies to the energy table for this edge.
///
/// @param infile - [in/out] - the binary input file, that has already been advanced
///   so that it's pointing at the position where the information for this
///   edge is stored.  At the end of this method, the infile will be advanced
///   past all the data for this edge.
///
void PDEdge::read_edge_energies_from_file( std::ifstream & infile )
{
std::cerr << "Reading Edge: " << get_node_index(0) << " " << get_node_index(1) << std::endl;
int node1_num_states_in_file = get_pd_node(0)->get_num_states_in_file();
int node2_num_states_in_file = get_pd_node(1)->get_num_states_in_file();

FArray1A_int node1_file_states_2_instance_states(
get_pd_node(0)->get_file_states_2_instance_states_array(),
node1_num_states_in_file);

FArray1A_int node2_file_states_2_instance_states(
get_pd_node(1)->get_file_states_2_instance_states_array(),
node2_num_states_in_file);

FArray1A_int node1_aatypes_file_states(
get_pd_node(0)->get_aatypes_for_file_states(),
node1_num_states_in_file);

FArray1A_int node2_aatypes_file_states(
get_pd_node(1)->get_aatypes_for_file_states(),
node2_num_states_in_file);

int num_file_aa = get_pdig_owner()->get_num_file_aatypes();
int num_aa = get_pdig_owner()->get_num_aatypes();

FArray1A_int node1_num_file_states_for_aa(
get_pd_node(0)->get_num_file_states_for_aa(),
num_file_aa);

FArray1A_int node2_num_file_states_for_aa(
get_pd_node(1)->get_num_file_states_for_aa(),
num_file_aa);

int sqr_file_aa = num_file_aa * num_file_aa;
FArray2D_bool aa_neighbors( num_aa, num_aa, false );
bool* aa_neighbor_buffer = new bool[ sqr_file_aa ];
for (int ii = 0; ii <= sqr_file_aa; ++ii)
aa_neighbor_buffer[ii] = false;

//std::cerr << "square file aa: " << sqr_file_aa << std::endl;
infile.read( (char*) aa_neighbor_buffer, sizeof( bool ) * sqr_file_aa );
//int num_bools_read = infile.gcount();
//debug_assert( num_bools_read == sizeof( bool ) * sqr_file_aa );

int buffer_index = 0;
int num_pair_energies = 0;
for (int ii = 1; ii <= num_file_aa; ++ii) {
for (int jj = 1; jj <= num_file_aa; ++jj) {
aa_neighbors(jj, ii) = aa_neighbor_buffer[ buffer_index ];
//std::cerr << aa_neighbor_buffer[ buffer_index ];
if ( aa_neighbor_buffer[ buffer_index ] ) {
num_pair_energies += node1_num_file_states_for_aa(ii) *
node2_num_file_states_for_aa( jj );
}
++buffer_index;

}
//std::cerr << std::endl;
}

//std::cerr << "num pair energies: " << num_pair_energies << std::endl;

two_body_energies_.set_sparse_aa_info( aa_neighbors );

float * energies_buffer = new float [ num_pair_energies ];

infile.read( (char*) energies_buffer, sizeof( float ) * num_pair_energies );
//int read_num_floats = infile.gcount();
buffer_index = 0;

for (int ii = 1; ii <= node1_num_states_in_file; ++ii) {
int ii_aa = node1_aatypes_file_states( ii );
int ii_instance_state = node1_file_states_2_instance_states( ii );
for (int jj = 1; jj <= node2_num_states_in_file; ++jj) {
int jj_aa = node2_aatypes_file_states( jj );
int jj_instance_state = node2_file_states_2_instance_states( jj );

if ( ! aa_neighbors(jj_aa, ii_aa ) ) continue;

debug_assert( buffer_index < num_pair_energies );
float energy = energies_buffer[ buffer_index ];
++buffer_index;

if ( ii_instance_state == -1 || jj_instance_state == -1 ) continue;

two_body_energies_.set(
get_pd_node(0)->get_sparse_mat_info_for_state(ii_instance_state),
get_pd_node(1)->get_sparse_mat_info_for_state(jj_instance_state),
energy);
}
}

delete [] aa_neighbor_buffer;
delete [] energies_buffer;
}

/// @brief Advances the infile past the section describing an edge between
/// one (or two) file nodes that do not correspond to an instance node.
///
/// This static method is very similar to the non-static
/// get_edge_energies_from_file() above, except that it requires
/// more input parameters, handed to it by a PDInteractionGraph instance
/// since the method cannot request information from two nodes of the graph.
/// instead of allocating a large buffer for the energies, it uses
/// the ofstream method seekg to advance the read head past the energies
/// for the absent edge.
///
/// @param infile - [in/out] - the binary input file, that's been advanced to point
///   at information for a file edge that will not be included in the instance
///   interaction graph.  At the end of this method, the file will point
///   just past the information for this edge.
/// @param num_file_aa - [in] - the number of aa types according to the file
/// @param node1_num_file_states_for_aa - [in] - num states for aa on node 1
/// @param node2_num_file_states_for_aa - [in] - num states for aa on node 2
///
void
PDEdge::skip_over_edge_energies_from_file
(
std::ifstream & infile,
int num_file_aa,
FArray1_int & node1_num_file_states_for_aa,
FArray1_int & node2_num_file_states_for_aa
)
{
int sqr_file_aa = num_file_aa * num_file_aa;
bool* aa_neighbor_buffer = new bool[ sqr_file_aa ];

infile.read( (char*) aa_neighbor_buffer, sizeof(bool) * sqr_file_aa);

int buffer_index = 0;
int num_pair_energies = 0;
for (int ii = 1; ii <= num_file_aa; ++ii) {
for (int jj = 1; jj <= num_file_aa; ++jj) {
if ( aa_neighbor_buffer[ buffer_index ] ) {
num_pair_energies += node1_num_file_states_for_aa( ii ) *
node2_num_file_states_for_aa( jj );
}
++buffer_index;
}
}

//skip forward
infile.seekg( sizeof( float ) * num_pair_energies, std::ios::cur );
delete [] aa_neighbor_buffer;
}

/// @brief Writes the energies for this edge to a binary output file.
///
/// First writes out the amino acid neighbor information for this edge
/// then writes out the energies.
///
/// @param outfile - [in/out] - the binary output file to write to.
///
void PDEdge::write_edge_energies_to_file( std::ofstream & outfile )
{
FArray1A_int node1_aatypes_for_state(
get_pd_node(0)->get_aatypes_for_states(),
get_num_states_for_node(0) );

FArray1A_int node2_aatypes_for_state(
get_pd_node(1)->get_aatypes_for_states(),
get_num_states_for_node(1) );

//write the aa-neighbor information
int num_aa = get_pdig_owner()->get_num_aatypes();
int sqr_num_aa = num_aa * num_aa;
bool * aa_neighbor_buffer = new bool[ sqr_num_aa ];
int buffer_index = 0;

for (int ii = 1; ii <= num_aa; ++ii) {
for (int jj = 1; jj <= num_aa; ++jj) {
if ( two_body_energies_.get_sparse_aa_info(ii, jj) ) {
aa_neighbor_buffer[ buffer_index ] = true;
} else {
aa_neighbor_buffer[ buffer_index ] = false;
}
//std::cerr << aa_neighbor_buffer[ buffer_index ];
++buffer_index;
}
//std::cerr << std::endl;
}
outfile.write( (char*) aa_neighbor_buffer, sizeof( bool ) * sqr_num_aa );

float * energy_buffer = new float[ two_body_energies_.get_table_size() ];
buffer_index = 0;
for (int ii = 1; ii <= get_num_states_for_node(0); ++ii ) {
int ii_aa = node1_aatypes_for_state( ii );
for (int jj = 1; jj <= get_num_states_for_node(1); ++jj ) {
int jj_aa = node2_aatypes_for_state( jj );
if ( ! two_body_energies_.get_sparse_aa_info(ii_aa, jj_aa) ) continue;

energy_buffer[ buffer_index ] = get_two_body_energy( ii, jj );
++buffer_index;
}
}

outfile.write( (char*) energy_buffer, sizeof(float) * two_body_energies_.get_table_size() );
std::cerr << "Writing edge: " << get_node_index(0) << " " << get_node_index(1);
std::cerr << "; num_energies: " << two_body_energies_.get_table_size() << std::endl;
delete [] aa_neighbor_buffer;
delete [] energy_buffer;
}

*/

/// @brief allow derived class to prep this class for simA, but guarantee no call to delete this;
void PDEdge::declare_energies_final_no_deletion()
{
	prepare_for_simulated_annealing_no_deletion();
}

/// @brief - allow derived class to prep this class for simA, but guarantee no call to delete this;
void PDEdge::prepare_for_simulated_annealing_no_deletion() //hook for derived classes
{
	if ( ! energies_updated_since_last_prep_for_simA_ ) return;

	//std::cout << "PDEdge: prepare_for_simulated_annealing_no_deletion(): two_body_energies table size before drop call: " << two_body_energies_.get_table_size() << std::endl;
	two_body_energies_.drop_zero_submatrices_where_possible();
	//std::cout << "PDEdge: prepare_for_simulated_annealing_no_deletion(): two_body_energies table size after drop call: " << two_body_energies_.get_table_size() << std::endl;
	energies_updated_since_last_prep_for_simA_ = false;
}

//// @brief - if edge table is empty, returns true -- assumes
///   prepare_for_simulated_annealing_no_deletion() has been called first.
bool PDEdge::pd_edge_table_all_zeros() const
{
	debug_assert( ! energies_updated_since_last_prep_for_simA_ );
	return ( two_body_energies_.get_table_size() == 0);
}


/// @brief drops any amino-acid neighbor submatrix of the two-body energy table
/// when the magnitudes of the energies stored in that submatrix do not exceed
/// the input parameter, epsilon.  Dropping submatrices that contain zero
/// energies is a special case of this function where epsilon == 0.
///
/// @param epsilon - [in] - the magnitude threshold for keeping amino-acid neighbor
/// submatrices.
///
void PDEdge::drop_small_submatrices_where_possible( float epsilon )
{
	two_body_energies_.drop_small_submatrices_where_possible( epsilon );
	return;
}

/// @brief drops amino-acid neighbor submatrices when they do not contain any non-zero
/// entries.  Represents a special case of drop_small_submatrices_where_possible
///
void PDEdge::drop_zero_submatrices_where_possible()
{
	drop_small_submatrices_where_possible( 0.0f );
	return;
}


/// @brief returns the memory usage of the two body energy table for this edge
///
int PDEdge::get_two_body_table_size() const
{
	return two_body_energies_.get_table_size();
}

//----------------------------------------------------------------------------//
//-------- Sparse Pairwise Decomposable Interaction Graph Class --------------//
//----------------------------------------------------------------------------//

/// @brief main constructor: no default nor copy constructors provided.
///
///
/// @param num_nodes - [in] - the number of nodes in this graph
PDInteractionGraph::PDInteractionGraph(int num_nodes) :
	PrecomputedPairEnergiesInteractionGraph( num_nodes ),
	num_aa_types_( -1 ), num_commits_since_last_update_(0),
	total_energy_current_state_assignment_(0),
	total_energy_alternate_state_assignment_(0),
	node_considering_alt_state_( -1 )
{}

void
PDInteractionGraph::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base )
{
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

	// determine max # of residue type groups
	Size max_nresgroups = 0;
	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size ii_nresgroups =  rot_sets.rotamer_set_for_moltenresidue( ii )->get_n_residue_groups();
		if ( ii_nresgroups > max_nresgroups ) max_nresgroups = ii_nresgroups;
	}

	//"aa types" means "distinct groups of rotamers" -- this ig has no idea
	// what an amino acid is or why they might be different from one another
	//set_num_aatypes( (int) max_nresgroups );
	num_aa_types_ = (int) max_nresgroups;

	for ( Size ii = 1; ii <= rot_sets.nmoltenres(); ++ii ) {
		Size const ii_num_states = rot_sets.rotamer_set_for_moltenresidue( ii )->num_rotamers();
		set_num_states_for_node( ii, ii_num_states );

		// figure out which residue-type group each rotamer is a member of
		std::vector< int > aatype_for_state( ii_num_states + 1, 0 );
		Size curr_resgroup = 1;
		Size count_for_resgroup = 1;
		Size const ii_nresgroups = rot_sets.rotamer_set_for_moltenresidue( ii )->get_n_residue_groups();
		for ( uint jj = 1; jj <= ii_num_states; ++jj ) {

			aatype_for_state[ jj ] = curr_resgroup;
			++count_for_resgroup;
			while ( count_for_resgroup > rot_sets.rotamer_set_for_moltenresidue( ii )->get_n_rotamers_for_residue_group( curr_resgroup ) ) {
				// increment curr_restype and skip over restypes with 0 rotamers
				++curr_resgroup;
				count_for_resgroup = 1;
				if ( curr_resgroup > ii_nresgroups ) break;
			}

		}

		get_pd_node( ii )->set_amino_acid_types( aatype_for_state );
	}

}


/// @details The actual meaning of the integers used to represent amino acid types
/// is not important, rather each state is labeled as being of some amino
/// acid type.  The amino acid types mean little more than equivalence classes;
/// two states are either in the same equivalence class or they are in
/// different ones.
/// this function should be called once and only once
///
/// @param num_aa_types - [in] - the number of amino acid types
/*void PDInteractionGraph::set_num_aatypes(int num_aa_types)
{
debug_assert( num_aa_types_ == -1 && num_aa_types > 0 );
num_aa_types_ = num_aa_types;
return;
}*/

/// @brief returns the number of different amino acid types
int  PDInteractionGraph::get_num_aatypes() const
{ return num_aa_types_;}

void
PDInteractionGraph::add_edge(int node1, int node2)
{
	InteractionGraphBase::add_edge( node1, node2 );
	PDEdge* newedge = get_pd_edge( node1, node2 );
	newedge->force_all_aa_neighbors();
}


/// @brief accessor
float
PDInteractionGraph::get_one_body_energy_for_node_state( int node, int state)
{
	return get_pd_node( node )->get_one_body_energy( state );
}

/// @brief assigns the state of all nodes in the interaction graph to their unassigned
/// or zero states.
void PDInteractionGraph::blanket_assign_state_0()
{
	//a state assignment of 0 means "unassigned".
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_pd_node(ii)->assign_zero_state();
	}
	total_energy_current_state_assignment_ = 0;
	return;
}

/// @brief sets the state on node node_ind to new_state
///
/// @param node_ind - [in] - the index of the node in question
/// @param new_state - [in] - the new state the node is being assigned to
float PDInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	get_pd_node( node_ind )->assign_state(new_state);
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
float PDInteractionGraph::set_network_state( ObjexxFCL::FArray1_int & node_states)
{
	//node_states.dimension( get_num_nodes() );
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_pd_node( ii )->assign_state( node_states(ii) );
	}
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

/// @brief considers altering the state of a particular node; returns the
/// change in energy that the state substitution would produce
///
/// To avoid too much numerical drift from accumulating, the bookkeeping
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
PDInteractionGraph::consider_substitution
(
	int node_ind,
	int new_state,
	float & delta_energy,
	float & prev_energy_for_node
)
{
	node_considering_alt_state_ = node_ind;
	delta_energy = get_pd_node( node_ind )->
		project_deltaE_for_substitution( new_state, prev_energy_for_node );

	//numerical drift accumulates in the following assignment
	total_energy_alternate_state_assignment_ =
		total_energy_current_state_assignment_ + delta_energy;

	return;
}

/// @brief Accepts (commits) the state change previously considered in a call to
/// consider_substitution and returns the energy of the entire graph
///
/// to avoid too much numerical drift from accumulating, the bookkeeping
/// arrays are updated once every 2^10 state commits
///
float
PDInteractionGraph::commit_considered_substitution()
{
	get_pd_node( node_considering_alt_state_ )->
		commit_considered_substitution();
	total_energy_current_state_assignment_ =
		total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals();
	}

	return total_energy_alternate_state_assignment_;
}

/// @brief removes all accumulated numerical drift and returns the
/// energy for the current state assignment.
float PDInteractionGraph::get_energy_current_state_assignment()
{
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}


float PDInteractionGraph::get_energy_PD_current_state_assignment()
{
	return total_energy_current_state_assignment_;
}

/// @details Iterates across nodes and then edges to look-up the energies
/// for the current state assignmnet removing any numerical drift which
/// accumulated in the member variable total_energy_current_state_assignment_.
void PDInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	//std::cerr << "updating internal energy totals: " << std::endl;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		//std::cerr << " ig_node " << ii << " = " << ((PDNode *) ig_nodes_[ ii ])
		// ->get_one_body_energy_current_state();

		total_energy_current_state_assignment_ += get_pd_node( ii )->
			get_one_body_energy_current_state();
	}

	//int counter = 0;
	for ( std::list<EdgeBase*>::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		//std::cerr << " ig_edge " << ++counter  << " =" <<
		//((PDEdge*) *iter)->get_current_two_body_energy();
		total_energy_current_state_assignment_ +=
			((PDEdge*) *iter)->get_current_two_body_energy();
	}

	//std::cerr << std::endl;

	num_commits_since_last_update_ = 0;
	return;
}

int PDInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for ( std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		sum += ((PDEdge*) *iter)->get_two_body_table_size();
	}
	return sum;
}

unsigned int
PDInteractionGraph::getMemoryUsageInBytes() const
{
	unsigned int total_memory = 0;
	total_memory += sizeof( PDInteractionGraph );
	//total_memory += InteractionGraphBase::getMemoryUsageInBytes();
	return total_memory;
}

void PDInteractionGraph::print_current_state_assignment() const
{
	std::cerr << "Curr States: ";
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		std::cerr << "(" << ii << ", ";
		std::cerr << get_pd_node(ii)->get_current_state() << ") ";
		get_pd_node(ii)->print_internal_energies();
	}
	std::cerr << std::endl;
}

/*
/// @brief allocates space for the arrays necessary for reading in from a file
///
void PDInteractionGraph::prepare_to_read_energies_from_file()
{
instance_node_2_file_node_.dimension( get_num_nodes() );
instance_node_2_file_node_ = -1;
}


/// @brief deallocates arrays used to read energies from a file
void PDInteractionGraph::declare_finished_reading_from_file()
{
instance_node_2_file_node_.dimension( 0 );
file_node_2_instance_node_.dimension( 0 );
aa_types_for_states_on_file_nodes_.dimension( 0 );


for (int ii = 1; ii <= get_num_nodes(); ++ii )
{
get_pd_node(ii)->clean_up_after_reading_energies_from_file();
}
}

/// @brief sets the number of amino acid types according to the input file
///
/// @param num_file_aatypes - [in] - the number of amino acid types
void PDInteractionGraph::set_num_file_aatypes( int num_file_aatypes )
{
num_file_aa_types_ = num_file_aatypes;
}

/// @brief returns the number of amino acid types according to the input file
///
int PDInteractionGraph::get_num_file_aatypes()
{
return num_file_aa_types_;
}

/// @brief sets the number of nodes for the graph described in the input file
///
/// @param num_nodes_in_file - [in] - the number of nodes in the file
void PDInteractionGraph::set_num_nodes_in_file( int num_nodes_in_file )
{
num_nodes_in_file_ = num_nodes_in_file;
file_node_2_instance_node_.dimension( num_nodes_in_file );
aa_types_for_states_on_file_nodes_.dimension( num_nodes_in_file );
num_file_states_for_aa_for_node_.dimension( num_nodes_in_file );

file_node_2_instance_node_ = -1;
}

/// @brief records the correspondence between an instance node and a file node
///
/// @param instance_node - [in] - the index of the instance node
/// @param file_node - [in] - the index of the file node
///
void PDInteractionGraph::set_node_correspondence
(
int instance_node,
int file_node
)
{
file_node_2_instance_node_( file_node ) = instance_node;
instance_node_2_file_node_( instance_node ) = file_node;
}


/// @brief records the number of states for a file node
///
/// If a file node corresponds to an instance node, then the instance
/// node keeps track of the information for the file node.  If the file
/// node doesn't correspond to any instance node, then the graph needs
/// to keep track of the file node's information.
/// Assumption: the correspondence between file nodes and instance nodes
/// has been completely specified before this method is called for the first time.
///
/// @param node - [in] - the index of the file node
/// @param num_file_states - [in] - the number of file states for that file node
///
void PDInteractionGraph::set_num_states_for_file_node
(
int node,
int num_file_states
)
{
if ( file_node_2_instance_node_( node ) == -1 ) {
aa_types_for_states_on_file_nodes_( node ).dimension( num_file_states );
num_file_states_for_aa_for_node_( node ).dimension( num_file_aa_types_ );
num_file_states_for_aa_for_node_( node ) = 0;
} else {
int instance_node = file_node_2_instance_node_( node );
get_pd_node(instance_node)->
prepare_to_read_energies_from_file( num_file_states );
}
}

/// @brief sets the amino acid type for a file state
///
/// @param node - [in] - the index of the file node
/// @param file_state - [in] - the file state
/// @param state_aa - [in] - the amino acid type for file_state
///
void PDInteractionGraph::set_aa_for_file_node_state
(
int node,
int file_state,
int state_aa
)
{
int instance_node = file_node_2_instance_node_( node );
if ( instance_node == -1 ) {
aa_types_for_states_on_file_nodes_( node )( file_state ) = state_aa;
++num_file_states_for_aa_for_node_( node )(state_aa);
} else {
get_pd_node( instance_node )->set_aa_for_file_state(file_state, state_aa);
}
}

/// @brief sets the correspondence between an instance state and a file state
///
/// @param node - [in] - the index of the instance node
/// @param state - [in] - the instance state
/// @param file_state - [in] - the file state that the instance state corresponds to
///
void PDInteractionGraph::set_correspondence_for_state(int node, int state, int file_state)
{
get_pd_node(node)->set_instance_state_correspondence( state, file_state );
}

/// @brief returns the file state that an instance state corresponds to
///
/// @param node - [in] - the instance node
/// @param state - [in] - the instance state
///
int PDInteractionGraph::get_correspondence_for_state(int node, int state )
{
return get_pd_node(node)->get_correspondence_for_state( state );
}

/// @brief returns true if an instance node corresponds to a node in the file
///
/// @param node - [in] - the instance node in question
///
bool PDInteractionGraph::get_node_corresponded_to_file_node( int node )
{
return get_pd_node(node)->get_node_corresponded_to_file_node( );
}

/// @brief returns the number of instance states not described by any file state for
/// for an instance node
///
/// @param node - [in] - the instance node in question
///
int PDInteractionGraph::get_num_rots_absent_from_file(int node)
{
return get_pd_node(node)->get_num_rots_absent_from_file();
}

/// @brief fills an input array with the indices of the instance states for an
/// instance node that did not correspond to any instance state of the input
/// file.
///
/// @param node - [in] - the instance node
/// @param rots_absent - [out] - the array in which to write the absent rotamers
///
void PDInteractionGraph::get_absent_rots(int node, FArray1_int & rots_absent )
{
get_pd_node(node)->get_absent_rots( rots_absent );
return;
}

/// @brief reads pair energies from a binary file
///
/// The edge energies in the binary input file have the following structure:
/// The number of edges in the graph is listed, and is followed by  a long list
/// of edges and edge information; each edge is identified by the indices
/// of the two file nodes that it is incident upon. Following the two
/// node indices, the rest of the information for that edge is stored.
/// The interaction graph reads in the file nodes and, if both nodes correspond
/// to instance nodes, creates a new edge between the instance nodes and
/// has the new edge read in the appropriate energies from the input file.
/// If either file node fails to correspond to some instance node, then the
/// graph invokes the static method PDEdge::skip_past_edge_energies() to advance
/// the read pointer past the information for this edge.
///
/// @param infile - [in/out] - the binary input file that has been advanced to point
///   at the number of edges in the file
///
void PDInteractionGraph::read_edge_energies_from_file( std::ifstream & infile )
{
int num_edges;
infile.read( (char*) &num_edges,4);

for (int ii = 1; ii <= num_edges; ++ii)
{
int node1, node2;
infile.read( (char*)  & node1, 4);
infile.read( (char*)   & node2, 4 );

int instance_node1 = file_node_2_instance_node_( node1 );
int instance_node2 = file_node_2_instance_node_( node2 );
if ( instance_node1 == -1 || instance_node2 == -1 )
{ //skip over the energies for this edge
if (instance_node1 == -1 ) {
FArray1A_int node1_num_states_for_aa(
num_file_states_for_aa_for_node_( node1 ), num_file_aa_types_);

if (instance_node2 == -1) {
FArray1A_int node2_num_states_for_aa(
num_file_states_for_aa_for_node_( node2 ), num_file_aa_types_);


PDEdge::skip_over_edge_energies_from_file( infile,
num_file_aa_types_,
node1_num_states_for_aa, node2_num_states_for_aa);
} else {
FArray1A_int node2_num_states_for_aa(
get_pd_node( instance_node2 )
->get_num_file_states_for_aa(), num_file_aa_types_ );

PDEdge::skip_over_edge_energies_from_file( infile,
num_file_aa_types_,
node1_num_states_for_aa, node2_num_states_for_aa);
}
} else {
FArray1A_int node1_num_states_for_aa(
get_pd_node( instance_node1 )
->get_num_file_states_for_aa(), num_file_aa_types_);

FArray1A_int node2_num_states_for_aa(
num_file_states_for_aa_for_node_( node2 ), num_file_aa_types_);

PDEdge::skip_over_edge_energies_from_file( infile,
num_file_aa_types_,
node1_num_states_for_aa, node2_num_states_for_aa);

}
} else {
add_edge( instance_node1, instance_node2 );
PDEdge* new_edge = (PDEdge*) find_edge( instance_node1, instance_node2);

new_edge->read_edge_energies_from_file( infile );
}

}
}

/// @brief writes the edge energies, in a binary format, to the outfile.
///
/// The graph writes out the number of edges, then proceeds to describe each
/// edge.  Each edge is identified by the two nodes its incident upon, and
/// each edge is responsible for outputting itself to the file.
///
/// @param outfile - [in/out] - the output file to write edge energies into
///
void PDInteractionGraph::write_edge_energies_to_file( std::ofstream & outfile )
{
for (int ii = 1; ii <= get_num_nodes(); ++ii) {
get_pd_node(ii)->prepare_to_write_to_file();
}

int num_edges = get_num_edges(); //O(N) size method
outfile.write( (char*) (&num_edges), 4);
for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
edge_iter != get_edge_list_end(); ++edge_iter) {
int first_node_ind = (*edge_iter)->get_first_node_ind();
int second_node_ind = (*edge_iter)->get_second_node_ind();
outfile.write( (char*) & first_node_ind, 4);
outfile.write( (char*) & second_node_ind, 4);
((PDEdge*) (*edge_iter) )->write_edge_energies_to_file( outfile );
}

for (int ii = 1; ii <= get_num_nodes(); ++ii) {
get_pd_node(ii)->clean_up_after_writing_to_file();
}

return;
}
*/

/// @details For instance in a graph with 6 vertices,{a,b,c,d,e,f}
/// a user may be interested in the sum of the one- and two-body energies
/// for vertices {a,b,c}.  The graph will return sum of the one body energies
/// for vertices a b and c and also any two-body energies for the edges in the
/// subgraph induced by a,b, and c.  (In this case, edges {a,b}, {a,c} and {b,c}
/// if these edges are part of the graph.  The edge {a,d} will not be counted
/// if it is part of the graph.).
/// ask the graph for the energies of the induced subgraph defined
/// by a particular group.
///
/// @param group_id - [in] - the groups for which you're interested in retrieving
///  energies of the induced subgraph
///
float
PDInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	float esum = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_pd_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id ) ) {
			esum += ((PDEdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}

	return esum;
}


unsigned int
PDInteractionGraph::count_static_memory() const {
	return sizeof( PDInteractionGraph );
}

unsigned int
PDInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = 0;
	total_memory += InteractionGraphBase::count_dynamic_memory();
	return total_memory;
}

/// @details PDInteractionGraph will return aa submatrices as requested.
bool
PDInteractionGraph::aa_submatrix_energies_retrievable() const
{
	return true;
}

int PDInteractionGraph::aatype_for_node_state(
	int node_ind,
	int node_state
) const
{
	return get_pd_node(node_ind )->aatype_for_state( node_state );
}

ObjexxFCL::FArray2D< core::PackerEnergy >
PDInteractionGraph::get_aa_submatrix_energies_for_edge(
	int node1,
	int node2,
	int node1aa,
	int node2aa
) const
{
	return get_pd_edge( node1, node2 )->get_aa_submatrix_energies( node1aa, node2aa );
}


/// @param node_index - [in] - the index of the node being created
/// @param num_states - [in] - the total number of states for the new node
///
NodeBase* PDInteractionGraph::create_new_node( int node_index, int num_states)
{
	PDNode* new_node = new PDNode(this, node_index, num_states);
	debug_assert( new_node != NULL );
	return new_node;
}


/// @param index1 - [in] - the smaller-indexed node this edge is incident upon
/// @param index2 - [in] - the larger-indexed node this edge is incident upon

EdgeBase* PDInteractionGraph::create_new_edge( int index1, int index2)
{
	return new PDEdge(this, index1, index2);
}


// <directed_design>


/// calls project_deltaE_for_substitution ( the weighted version )


void PDNode::project_deltaE_for_substitution
(
	int alternate_state,
	float & deltaE_unweighted,
	float & prevE_unweighted,
	float & deltaE_weighted,
	float & prevE_weighted,
	ObjexxFCL::FArray2D< core::PackerEnergy > const & weights
)
{

	// Step 1 - calculate unweighted energy, with side effect of filling alternate_state_one/two_body_energy

	deltaE_unweighted = project_deltaE_for_substitution(alternate_state,prevE_unweighted); // !!! need to return this parameter since the graph is adding it to total

	// Step 2 - calculate weighted energy on basis of info in:
	// - curr_state_one_body_energy_
	// - curr_state_two_body_energies_
	// - alternate_state_one_body_energy_
	// - alternate_state_two_body_energies_

	float curr_total_energy_weighted      = curr_state_one_body_energy_;
	float alternate_total_energy_weighted = alternate_state_one_body_energy_;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {

		//   const int i = get_node_index();
		//   const int j = get_adjacent_node(ii)->get_node_index();
		//const float bias_ii = get_bias(bias,i,j);//get_node_index(),get_adjacent_node(ii)->get_node_index());

		const float bias_ii = weights(get_incident_edge(ii)->get_first_node_ind(), // !!! symmetric bias matrix
			get_incident_edge(ii)->get_second_node_ind());

		curr_total_energy_weighted      += bias_ii * curr_state_two_body_energies_[ii];
		alternate_total_energy_weighted += bias_ii * alternate_state_two_body_energies_[ii];

	} // ii

	prevE_weighted = curr_total_energy_weighted;
	deltaE_weighted = alternate_total_energy_weighted - curr_total_energy_weighted;

	//std::cout << __func__ << "deltaE_unweighted = " << deltaE_unweighted << "\tdeltaE_weighted = " << deltaE_weighted << std::endl;

} // project_deltaE_for_substitution

// THIS IS DUPLCATED CODE: EdgeBase already provides this functionality
namespace {

int get_other_index(const EdgeBase* edge_base, const int index) {
	if ( index == edge_base->get_first_node_ind() ) {
		return edge_base->get_second_node_ind();
	} else if ( index == edge_base->get_second_node_ind() ) {
		return edge_base->get_first_node_ind();
	} else {
		TR.Fatal << "Got an index of " << index << " expected either " << edge_base->get_first_node_ind() << " or " << edge_base->get_second_node_ind() << std::endl;
		utility_exit_with_message("get_other_index(const EdgeBase* edge_base, const int index)  failed");
		return -1;
	}
} // get_other_index

} // namespace

/// @brief Paul's code

float
PDNode::get_weighted_energy_with_higher_indexed_nodes( ObjexxFCL::FArray2D< core::PackerEnergy > const & weights ) const
{
	float rval = curr_state_one_body_energy_;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		if ( get_other_index(get_incident_edge(ii),get_node_index()) > ii ) {
			const float weight_ii = weights(get_incident_edge(ii)->get_first_node_ind(), // !!! symmetric bias matrix
				get_incident_edge(ii)->get_second_node_ind());
			rval += weight_ii * curr_state_two_body_energies_[ii];
		}
	} // ii

	return rval;
} // PDNode::get_biased_energy_with_higher_indexed_nodes


float PDInteractionGraph::get_weighted_energy(const ObjexxFCL::FArray2D< core::PackerEnergy >& weights) const {
	// Compute and return biased energy
	float rval = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		rval += get_pd_node(ii)->get_weighted_energy_with_higher_indexed_nodes(weights);
	}
	return rval;
} // PDInteractionGraph::get_weighted_energy


float PDInteractionGraph::set_network_state( ObjexxFCL::FArray1_int & node_states, ObjexxFCL::FArray2D< core::PackerEnergy > const & weights) {
	set_network_state(node_states);
	return get_weighted_energy(weights);
}

float PDInteractionGraph::commit_considered_substitution(const ObjexxFCL::FArray2D< core::PackerEnergy >& weights) {
	commit_considered_substitution();
	return get_weighted_energy(weights);
}

void PDInteractionGraph::consider_substitution
(
	int node_ind,
	int new_state,
	float & deltaE_unweighted,
	float & prevE_unweighted,
	float & deltaE_weighted,
	float & prevE_weighted,
	const ObjexxFCL::FArray2D< core::PackerEnergy >& weights
)
{

	node_considering_alt_state_ = node_ind;
	get_pd_node( node_ind )->project_deltaE_for_substitution( new_state, deltaE_unweighted, prevE_unweighted, deltaE_weighted, prevE_weighted, weights );

	//numerical drift accumulates in the following assignment
	// !!! need to make sure to add the unweighted energy to this variable
	total_energy_alternate_state_assignment_ =
		total_energy_current_state_assignment_ + deltaE_unweighted;

	return;
}

// </directed_design>

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core
