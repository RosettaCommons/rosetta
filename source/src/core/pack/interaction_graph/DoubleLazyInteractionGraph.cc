// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/DoubleLazyInteractionGraph.cc
/// @brief  Interaction graph that computes each rotamer pair energy at most once
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <core/pack/interaction_graph/DoubleLazyInteractionGraph.hh>

// Utility headers
#include <utility/in_place_list.hh>
#include <utility/pointer/owning_ptr.hh>

/// C++ headers
#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// main constructor, no default or copy constructors
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
DoubleLazyNode::DoubleLazyNode(
	InteractionGraphBase * owner,
	int node_id,
	int num_states
) :
	OnTheFlyNode( owner, node_id, num_states ),
	current_state_( 0 ),
	curr_state_one_body_energy_( 0.0f ),
	curr_state_total_energy_( 0.0f ),
	alternate_state_( 0 ),
	alternate_state_one_body_energy_( 0 ),
	alternate_state_total_energy_( 0 ),
	alternate_state_is_being_considered_( false ),
	procrastinated_( false )
{
}

DoubleLazyNode::~DoubleLazyNode()
{}

void
DoubleLazyNode::prepare_for_simulated_annealing()
{
	if (! get_edge_vector_up_to_date() ) update_internal_vectors();
	/*for (int ii = 1; ii <= get_num_states(); ++ii) {
		mark_coordinates_current( ii ); /// What did this used to to?
	}*/
	return;
}

void
DoubleLazyNode::print() const
{

	std::cout << "DoubleLazyNode " << get_node_index() << " with " << get_num_states() << " states" << std::endl;
	std::cout << "curr_state " << current_state_ << " ";
	std::cout << "curr_state_sparse_mat_info_ ";
	std::cout << curr_state_sparse_mat_info_.get_aa_type() << " ";
	std::cout << curr_state_sparse_mat_info_.get_state_ind_for_this_aa_type() << " ";
	std::cout << "Curr One Body Energy: " << curr_state_one_body_energy_ << std::endl;
	std::cout << "Curr Two Body Energies:";
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii) {
		std::cout << " " << get_index_of_adjacent_node(ii) << ":" << curr_state_two_body_energies_[ ii ];
	}
	std::cout << std::endl;

	if ( ! alternate_state_is_being_considered_ ) return;
	std::cout << "Alt One Body Energy: " << alternate_state_one_body_energy_ << std::endl;
	std::cout << "Alt Two Body Energies:";
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii) {
		std::cout << " " << get_index_of_adjacent_node(ii) << ":" << alternate_state_two_body_energies_[ ii ];
	}
	std::cout << std::endl  << "-----------------" << std::endl;


}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////

unsigned int
DoubleLazyNode::count_static_memory() const
{
	return sizeof( DoubleLazyNode );
}

unsigned int
DoubleLazyNode::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyNode::count_dynamic_memory();

	total_memory += neighbors_curr_state_.size() * sizeof (int );

	//total_memory += aa_offsets_for_edges_.size() * sizeof( int );
	//total_memory += num_states_for_aa_type_for_higher_indexed_neighbor_.size() * sizeof( int );
	//total_memory += edge_matrix_ptrs_.size() * sizeof( FArray1Da< core::PackerEnergy > );
	total_memory += neighbors_curr_state_.size() * sizeof( int );
	total_memory += neighbors_curr_state_sparse_info_.size() * sizeof( SparseMatrixIndex );

	total_memory += curr_state_two_body_energies_.size() * sizeof( core::PackerEnergy );
	total_memory += alternate_state_two_body_energies_.size() * sizeof( core::PackerEnergy );

	return total_memory;
}

int
DoubleLazyNode::aatype_for_state( int state ) const
{
	return get_sparse_mat_info_for_state( state ).get_aa_type();
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::assign_zero_state()
{
	current_state_ = 0;
	alternate_state_ = 0;
	alternate_state_is_being_considered_ = false;

	curr_state_one_body_energy_ = core::PackerEnergy( 0.0 );
	//fills from [1] to end
	std::vector< core::PackerEnergy >::iterator position1 = curr_state_two_body_energies_.begin();
	++position1;
	std::fill( position1, curr_state_two_body_energies_.end(), core::PackerEnergy( 0.0 ));
	curr_state_total_energy_ = core::PackerEnergy( 0.0 );

	for (int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_dlazy_edge(ii)->acknowledge_state_zeroed( get_node_index() );
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::assign_state(int new_state)
{
debug_assert( new_state >= 0 && new_state <= get_num_states());

	if (new_state == 0) {
		assign_zero_state();
	} else {
		//std::cout << "assign_state: node -  " << get_node_index() <<
		// " new state " << new_state << "...";
		current_state_ = new_state;
		curr_state_sparse_mat_info_ = get_sparse_mat_info_for_state( current_state_ );
		curr_state_one_body_energy_ = get_one_body_energy( current_state_ );
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for (int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_dlazy_edge(ii)->acknowledge_state_change(
				get_node_index(),
				current_state_,
				curr_state_sparse_mat_info_,
				curr_state_two_body_energies_[ii]);

			curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
		}
		//std::cout<< "..done" << std::endl;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::partial_assign_state( int new_state )
{
	if (new_state == 0 ) {
		assign_zero_state();
		return;
	}

	current_state_ = new_state;
	curr_state_sparse_mat_info_ =
		get_sparse_mat_info_for_state( current_state_ );
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_dlazy_edge(ii)->acknowledge_partial_state_change(
			get_node_index(),
			current_state_,
			curr_state_sparse_mat_info_);
	}
	alternate_state_is_being_considered_ = false;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void DoubleLazyNode::complete_state_assignment()
{
	if ( current_state_ == 0 ) return;

	curr_state_total_energy_ = curr_state_one_body_energy_ =
		get_one_body_energy( current_state_ );
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii) {
		curr_state_two_body_energies_[ ii ] =
			get_incident_dlazy_edge( ii )->
			get_energy_following_partial_state_assignment();
		curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	}
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::commit_considered_substitution()
{
debug_assert( alternate_state_is_being_considered_ );

	current_state_ = alternate_state_;
	curr_state_sparse_mat_info_ = alt_state_sparse_mat_info_;
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

	//if ( procrastinated_ )
	//{
	//	curr_state_total_energy_ = curr_state_one_body_energy_;
	//	for (int ii = 1; ii <= get_num_incident_edges(); ++ii)
	//	{
	//		if ( curr_state_two_body_energies_[ ii ] == DoubleLazyEdge::NOT_YET_COMPUTED_ENERGY )
	//		{
	//			curr_state_two_body_energies_[ ii ] = compute_pair_energy_for_current_state( ii );
	//		}
	//		curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	//	}
	//	procrastinated_ = false;
	//	++num_procrastinated_committed;
	//}


	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_dlazy_edge(ii)->acknowledge_substitution(
			get_node_index(),
			alternate_state_two_body_energies_[ii],
			current_state_,
			curr_state_sparse_mat_info_
		);
	}

	alternate_state_is_being_considered_ = false;
	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyNode::compute_pair_energy_for_current_state(
	int edge_making_energy_request
)
{

	return compute_rotamer_pair_energy(
		edge_making_energy_request,
		current_state_,
		neighbors_curr_state_[ edge_making_energy_request ]
	);
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::acknowledge_neighbors_partial_state_substitution(
	int edge_to_altered_neighbor,
	int other_node_new_state,
	SparseMatrixIndex const & other_node_new_state_sparse_info
)
{
	curr_state_total_energy_ = 0;
	curr_state_two_body_energies_[ edge_to_altered_neighbor ] = 0;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	neighbors_curr_state_sparse_info_[ edge_to_altered_neighbor ]  =
		other_node_new_state_sparse_info;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::print_internal_energies() const
{
	std::cout << "curr_state " << current_state_ << " ";
	std::cout << "curr_state_sparse_mat_info_ ";
	std::cout << curr_state_sparse_mat_info_.get_aa_type() << " ";
	std::cout << curr_state_sparse_mat_info_.get_state_ind_for_this_aa_type() << " ";
	std::cout << "curr_state_one_body_energy_ ";
	std::cout << curr_state_one_body_energy_ << " ";
	std::cout << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii) {
		std::cout << "(" << curr_state_two_body_energies_[ ii ] << ") ";
	}
	std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyNode::update_internal_energy_sums()
{
debug_assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = 0;
	for (int ii = 1; ii <= get_num_incident_edges(); ++ii) {
		curr_state_total_energy_ += get_incident_dlazy_edge(ii)->get_current_two_body_energy();
	}
	curr_state_total_energy_ += curr_state_one_body_energy_;
	return;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void DoubleLazyNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	neighbors_curr_state_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_sparse_info_.resize( get_num_incident_edges() + 1);

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		neighbors_curr_state_sparse_info_[ii].set_aa_type( 1 );
	}

	curr_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	alternate_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	return;
}

// @ brief - allow derived class to "drive" through the deltaE calculation
void
DoubleLazyNode::calc_deltaEpd( int alternate_state )
{
	core::PackerEnergy dummy(0.0f);
	project_deltaE_for_substitution( alternate_state, dummy );
}


//-----------------------------------------------------------------//

core::PackerEnergy const DoubleLazyEdge::NOT_YET_COMPUTED_ENERGY = -1234;

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
DoubleLazyEdge::DoubleLazyEdge(
	InteractionGraphBase * owner,
	int first_node_ind,
	int second_node_ind
):
	OnTheFlyEdge( owner, first_node_ind, second_node_ind),
	sparse_aa_neighbors_(
		get_dlazy_ig_owner()->get_num_aatypes(),
		get_dlazy_ig_owner()->get_num_aatypes(),
		(unsigned char) 0 ),
	two_body_energies_(
		get_dlazy_ig_owner()->get_num_aatypes(),
		get_dlazy_ig_owner()->get_num_aatypes(),
		( ObjexxFCL::FArray2D< core::PackerEnergy > * ) ( 0 ) ),
	curr_state_energy_( 0.0f ),
	partial_state_assignment_( false ),
	ran_annealing_since_pair_energy_table_cleared_( false ),
	// by default, do not use the memory-capping behavior of the DLIG;
	// an index >= 0 means use memory-capping behavior
	edge_index_( -1 )
{
}

/// @details explicity deallocate the FArrays that were dynamically
/// allocated
DoubleLazyEdge::~DoubleLazyEdge()
{
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			delete two_body_energies_( jj, ii );
		}
	}
}

void
DoubleLazyEdge::set_sparse_aa_info(
	ObjexxFCL::FArray2_bool const & aa_neighbors
)
{
	sparse_aa_neighbors_ = aa_neighbors;
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			if ( two_body_energies_( jj, ii ) && sparse_aa_neighbors_( jj, ii )) {
				(*two_body_energies_(jj,ii)) = NOT_YET_COMPUTED_ENERGY;
			} else if ( two_body_energies_(jj,ii) ) {
				delete two_body_energies_(jj,ii); two_body_energies_(jj,ii) = 0;
			}
		}
	}
	ran_annealing_since_pair_energy_table_cleared_ = false;
}

/// @brief returns whether two amino acid types are represented as neighbors
bool DoubleLazyEdge::get_sparse_aa_info( int node1aa, int node2aa ) const
{
	return sparse_aa_neighbors_( node2aa, node1aa );
}

void DoubleLazyEdge::force_aa_neighbors(int node1aa, int node2aa)
{
	sparse_aa_neighbors_( node2aa, node1aa ) = 1;
}

void DoubleLazyEdge::force_all_aa_neighbors()
{
	sparse_aa_neighbors_ = 1;
	ran_annealing_since_pair_energy_table_cleared_ = false;
}

core::PackerEnergy DoubleLazyEdge::get_two_body_energy(
	int const node1state,
	int const node2state
) const
{
	SparseMatrixIndex node1info = get_otf_node(0)->get_sparse_mat_info_for_state( node1state );
	SparseMatrixIndex node2info = get_otf_node(1)->get_sparse_mat_info_for_state( node2state );
	return get_two_body_energy_smi( node1state, node2state, node1info, node2info );
}

core::PackerEnergy
DoubleLazyEdge::get_two_body_energy_smi(
	int const node1state,
	int const node2state,
	SparseMatrixIndex const & node1info,
	SparseMatrixIndex const & node2info
) const
{
	if ( ! sparse_aa_neighbors_( node2info.get_aa_type(), node1info.get_aa_type() )) {
		return 0.0;
	}
	prep_aa_submatrix( node1info.get_aa_type(), node2info.get_aa_type() );
	core::PackerEnergy energy = read_aa_submatrix( node1info, node2info );

	if ( energy == NOT_YET_COMPUTED_ENERGY ) {
		energy = get_otf_node(0)->compute_rotamer_pair_energy(
			get_edges_position_in_nodes_edge_vector( 0 ),
			node1state, node2state ) * edge_weight();
		set_aa_submatrix( node1info, node2info, energy );
	}
	return energy;
}

/// @details -- const to be called within const functions; updates mutable data
void
DoubleLazyEdge::prep_aa_submatrix(
	int node1aa,
	int node2aa
) const
{
debug_assert( sparse_aa_neighbors_( node2aa, node1aa ));

	if ( ! two_body_energies_( node2aa, node1aa ) ) {
		two_body_energies_( node2aa, node1aa ) = new ObjexxFCL::FArray2D< core::PackerEnergy >(
			get_dlazy_node(1)->get_num_states_for_aa_types()[ node2aa ],
			get_dlazy_node(0)->get_num_states_for_aa_types()[ node1aa ],
			NOT_YET_COMPUTED_ENERGY );
		if ( edge_index_ != -1 ) {
			get_dlazy_ig_owner()->note_submatrix_added(
				edge_index_, submatrix_index( node1aa, node2aa ), submatrix_size( node1aa, node2aa ));
		}
	}
}

core::PackerEnergy
DoubleLazyEdge::read_aa_submatrix(
	SparseMatrixIndex node1info,
	SparseMatrixIndex node2info
) const
{
debug_assert( two_body_energies_( node2info.get_aa_type(), node1info.get_aa_type() ) );

	return (*two_body_energies_( node2info.get_aa_type(), node1info.get_aa_type() ))
		( node2info.get_state_ind_for_this_aa_type(), node1info.get_state_ind_for_this_aa_type() );
}

/// @details updates mutable data
void
DoubleLazyEdge::set_aa_submatrix(
	SparseMatrixIndex node1info,
	SparseMatrixIndex node2info,
	core::PackerEnergy setting
) const
{
debug_assert( two_body_energies_( node2info.get_aa_type(), node1info.get_aa_type() ) );
	(*two_body_energies_( node2info.get_aa_type(), node1info.get_aa_type() ))
		( node2info.get_state_ind_for_this_aa_type(), node1info.get_state_ind_for_this_aa_type() ) = setting;
}

void
DoubleLazyEdge::declare_energies_final()
{}

void
DoubleLazyEdge::prepare_for_simulated_annealing()
{
	if ( ! ran_annealing_since_pair_energy_table_cleared_ ) {
		ran_annealing_since_pair_energy_table_cleared_ = true;
		for ( Size ii = 1; ii <= sparse_aa_neighbors_.size2(); ++ii ) {
			for ( Size jj = 1; jj <= sparse_aa_neighbors_.size1(); ++jj ) {
				if ( sparse_aa_neighbors_( jj, ii ) ) {
					return;
				}
			}
		}
		delete this;
	}
}

unsigned int
DoubleLazyEdge::count_static_memory() const
{
	return sizeof( DoubleLazyEdge );
}

/// @details add up all the submatrices
unsigned int
DoubleLazyEdge::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyEdge::count_dynamic_memory();
	total_memory += sparse_aa_neighbors_.size() * sizeof( unsigned char );
	total_memory += actual_twobody_memory_use();

	/// account for the memory use of all the pointers and the FArrays themselves (not their data)
	total_memory += two_body_energies_.size() * sizeof( ObjexxFCL::FArray2D< core::PackerEnergy > * );
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			if ( two_body_energies_(jj,ii) ) {
				total_memory += sizeof( ObjexxFCL::FArray2D< core::PackerEnergy > );
			}
		}
	}
	return total_memory;
}

unsigned int
DoubleLazyEdge::actual_twobody_memory_use() const
{
	unsigned int total_memory( 0 );
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			if ( two_body_energies_(jj,ii) ) {
				total_memory += two_body_energies_(jj,ii)->size() * sizeof( core::PackerEnergy );
			}
		}
	}
	return total_memory;
}


unsigned int
DoubleLazyEdge::potential_twobody_memory_use() const
{
	unsigned int total_memory( 0 );
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			if ( sparse_aa_neighbors_(jj,ii) ) {
				total_memory +=
					get_dlazy_node(1)->get_num_states_for_aa_types()[ jj ] *
					get_dlazy_node(0)->get_num_states_for_aa_types()[ ii ] *
					sizeof( core::PackerEnergy );
			}
		}
	}
	return total_memory;
}

/// @details This function is only invoked if the graph is going to attempt to use
/// a memory-use cap.  The edge_index_ serves to both identify that the memory-use
/// cap is in effect, and to allow efficient communication between the graph and
/// its edges.
void
DoubleLazyEdge::set_edge_index(
	int index
)
{
debug_assert( index > -1 && edge_index_ == -1 ); // set this only once
	edge_index_ = index;
}

int
DoubleLazyEdge::drop_aa_submatrix(
	int submat_ind
) const
{
	std::pair< int, int > aainds = aa_indices_from_submatrix_index( submat_ind );
debug_assert( two_body_energies_( aainds.second, aainds.first ) );
	delete two_body_energies_( aainds.second, aainds.first ); two_body_energies_( aainds.second, aainds.first ) = 0;
	return submatrix_size( aainds.first, aainds.second );
}


void
DoubleLazyEdge::set_edge_weight( Real weight )
{
	Real const scale_factor = weight / edge_weight();
	for ( Size ii = 1; ii <= two_body_energies_.size2(); ++ii ) {
		for ( Size jj = 1; jj <= two_body_energies_.size1(); ++jj ) {
			if ( two_body_energies_(jj,ii) ) {
				ObjexxFCL::FArray2D< core::PackerEnergy > & ii_jj_table( *two_body_energies_(jj,ii)  );
				for ( Size kk = 0; kk < ii_jj_table.size(); ++kk ) {
					if ( ii_jj_table[ kk ] != NOT_YET_COMPUTED_ENERGY ) {
						ii_jj_table[ kk ] *= scale_factor;
					}
				}
			}
		}
	}
	edge_weight( weight );
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyEdge::get_current_two_body_energy() const
{
	return curr_state_energy_;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyEdge::acknowledge_state_change(
	int node_ind,
	int new_state,
	SparseMatrixIndex const & new_state_sparse_info,
	core::PackerEnergy & new_energy
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	int nodes_curr_states[2];
	SparseMatrixIndex nodes_curr_states_sparse_info[2];

	nodes_curr_states[ node_substituted ] = new_state;
	nodes_curr_states_sparse_info[ node_substituted ] = new_state_sparse_info;

	nodes_curr_states[ node_not_substituted ] =
		get_dlazy_node( node_not_substituted )->get_current_state();
	nodes_curr_states_sparse_info[ node_not_substituted ] =
		get_dlazy_node( node_not_substituted )->
		get_sparse_mat_info_for_curr_state();

	get_energy_for_state_pair( nodes_curr_states, nodes_curr_states_sparse_info );
	new_energy = curr_state_energy_;

	get_dlazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		new_state,
		new_state_sparse_info
	);

	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyEdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;
	SparseMatrixIndex dummy_sparse_info;
	dummy_sparse_info.set_aa_type( 1 );
	dummy_sparse_info.set_state_ind_for_this_aa_type(1);

	get_dlazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		0,
		dummy_sparse_info
	);
	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void DoubleLazyEdge::acknowledge_partial_state_change(
	int node_ind,
	int new_state,
	SparseMatrixIndex const & new_state_sparse_info
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;

	get_dlazy_node( node_not_substituted )->acknowledge_neighbors_partial_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		new_state,
		new_state_sparse_info
	);
	partial_state_assignment_ = true;
	return;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyEdge::get_energy_following_partial_state_assignment()
{
	if (partial_state_assignment_) {
		int nodes_curr_states[2];
		SparseMatrixIndex nodes_curr_states_sparse_info[2];

		for (int ii = 0; ii < 2; ++ii) {
			nodes_curr_states[ ii ] =
				get_dlazy_node( ii )->get_current_state();
			nodes_curr_states_sparse_info[ ii ] =
				get_dlazy_node( ii )->get_sparse_mat_info_for_curr_state();
		}
		get_energy_for_state_pair( nodes_curr_states, nodes_curr_states_sparse_info );
		partial_state_assignment_ = false;
	}
	return curr_state_energy_;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
int DoubleLazyEdge::get_two_body_table_size() const
{
	return 0; // TEMP!!! two_body_energies_.get_table_size();
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyEdge::get_energy_for_state_pair(
	int nodes_states[ 2 ],
	SparseMatrixIndex sparse_matrix_indices[ 2 ]
)
{
	bool one_node_in_zero_state = ( nodes_states[0] == 0 || nodes_states[1] == 0 );

	if (  one_node_in_zero_state ) {
		curr_state_energy_ = 0;
	} else {
		if ( ! sparse_aa_neighbors_( sparse_matrix_indices[1].get_aa_type(), sparse_matrix_indices[0].get_aa_type() ) ) {
			curr_state_energy_ = 0;
		} else {
			curr_state_energy_ = get_two_body_energy_smi(
				nodes_states[0], nodes_states[1],
				sparse_matrix_indices[0], sparse_matrix_indices[1] );
		}
	}
}

ObjexxFCL::FArray2D< core::PackerEnergy >
DoubleLazyEdge::get_aa_submatrix_energies(
	int node1aa,
	int node2aa
) const
{
	if ( ! sparse_aa_neighbors_( node2aa, node1aa )) {
		ObjexxFCL::FArray2D< core::PackerEnergy > empty;
		return empty;
	}

	prep_aa_submatrix( node1aa, node2aa );

debug_assert( two_body_energies_( node2aa, node1aa ) );
	ObjexxFCL::FArray2D< core::PackerEnergy > submat( * two_body_energies_( node2aa, node1aa ));
	int const iioffset = get_dlazy_node(0)->get_state_offset_for_aatype( node1aa );
	int const jjoffset = get_dlazy_node(1)->get_state_offset_for_aatype( node2aa );

	for ( Size ii = 1; ii <= submat.size2(); ++ii ) {
		for ( Size jj = 1; jj <= submat.size1(); ++jj ) {
			if ( submat(jj,ii) == NOT_YET_COMPUTED_ENERGY ) {
				int const ii_state = ii + iioffset;
				int const jj_state = jj + jjoffset;
				core::PackerEnergy iijjenergy = get_otf_node(0)->compute_rotamer_pair_energy(
					get_edges_position_in_nodes_edge_vector( 0 ),
					ii_state, jj_state );
				(*two_body_energies_( node2aa, node1aa ))( jj, ii ) = iijjenergy;
				submat(jj,ii) = iijjenergy;
			}
		}
	}

	if ( edge_index_ != -1 ) {
		get_dlazy_ig_owner()->note_submatrix_accessed( edge_index_, submatrix_index( node1aa, node2aa ) );
	}

	return submat;
}


/*void
DoubleLazyEdge::wipe_two_body_energies_for_node_state(
	int node,
	int state
)
{
	int other_node = node == 0 ? 1 : 0;

	SparseMatrixIndex states[ 2 ];
	states[ node ] = get_dlazy_node( node )->get_sparse_mat_info_for_state( state );
	int const other_node_num_states = get_dlazy_node( other_node )->get_num_states();

	for (int ii = 1; ii <= other_node_num_states; ++ii) {
		states[ other_node ] = get_dlazy_node( other_node )->get_sparse_mat_info_for_state( ii );
		two_body_energies_.set( states[ 0 ], states[ 1 ], NOT_YET_COMPUTED_ENERGY );

	}

}*/

void
DoubleLazyEdge::print_current_energy() const
{
	std::cout << "DoubleLazyEdge: " << get_node_index( 0 ) << "/" << get_node_index( 1 );
	std::cout << " energy= " << curr_state_energy_ << std::endl;
}

/// @details -- node1aa and node2aa are both indexed from 1
int
DoubleLazyEdge::submatrix_index( int node1aa, int node2aa ) const
{
	return ( node1aa - 1 ) * get_dlazy_ig_owner()->get_num_aatypes() + node2aa - 1;
}

int
DoubleLazyEdge::submatrix_size( int node1aa, int node2aa ) const
{
	return get_dlazy_node(1)->get_num_states_for_aa_types()[ node2aa ] *
		get_dlazy_node(0)->get_num_states_for_aa_types()[ node1aa ] * sizeof( core::PackerEnergy );
}

std::pair< int, int >
DoubleLazyEdge::aa_indices_from_submatrix_index( int submat_ind ) const
{
	int node1aa = submat_ind / get_dlazy_ig_owner()->get_num_aatypes() + 1;
	int node2aa = submat_ind % get_dlazy_ig_owner()->get_num_aatypes() + 1;
	return std::make_pair( node1aa, node2aa );

}


//-------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
DoubleLazyInteractionGraph::DoubleLazyInteractionGraph(
	int numNodes
) :
	OnTheFlyInteractionGraph( numNodes ),
	memory_max_for_rpes_( 0 ),
	curr_memory_for_rpes_( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
DoubleLazyInteractionGraph::~DoubleLazyInteractionGraph()
{
#ifdef APL_FULL_DEBUG
	unsigned int edge_potential_usage( 0 );
	unsigned int edge_actual_usage( 0 );
	for( std::list< EdgeBase* >::const_iterator eiter = get_edge_list_begin(),
			eiter_end = get_edge_list_end(); eiter != eiter_end; ++eiter ) {
		DoubleLazyEdge const * dledge = static_cast< DoubleLazyEdge const * >  (* eiter );
		edge_potential_usage += dledge->potential_twobody_memory_use();
		edge_actual_usage += dledge->actual_twobody_memory_use();
	}
	std::cout << "~DoubleLazyIG: used " << edge_actual_usage << " bytes instead of " << edge_potential_usage << std::endl;
#endif
}

void
DoubleLazyInteractionGraph::initialize(
	rotamer_set::RotamerSetsBase const & rot_sets_base
)
{
	parent::initialize( rot_sets_base );
	sqr_num_aa_types_ = get_num_aatypes() * get_num_aatypes();
}

void
DoubleLazyInteractionGraph::prepare_for_simulated_annealing()
{
	InteractionGraphBase::prepare_for_simulated_annealing();

	if ( memory_max_for_rpes_ == 0 ) return;

	/// Ok -- now allocate space for our in-place list of edge submatrix-access history information
	int const num_edges = get_num_edges();
	int const n_submatrices =  num_edges * sqr_num_aa_types_;
	dlazy_edge_vector_.resize( num_edges );
	aa_submatrix_history_list_ = InPlaceIntListOP( new utility::in_place_list< int >( n_submatrices ) );

	int count = 0; // index nodes from 0 -- makes module arithmetic easier
	for( std::list< EdgeBase* >::const_iterator eiter = get_edge_list_begin(),
			eiter_end = get_edge_list_end(); eiter != eiter_end; ++eiter ) {
		DoubleLazyEdge * dledge = static_cast< DoubleLazyEdge * >  (* eiter );
		dledge->set_edge_index( count );
		dlazy_edge_vector_[ count ] = dledge;
		++count;
	}
}


/// @details
void
DoubleLazyInteractionGraph::blanket_assign_state_0()
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_dlazy_node( ii )->assign_zero_state();
	}
	total_energy_current_state_assignment_ = 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	get_dlazy_node( node_ind )->assign_state( new_state );
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyInteractionGraph::set_network_state( ObjexxFCL::FArray1_int & node_states)
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_dlazy_node( ii )->partial_assign_state( node_states( ii ) );
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_dlazy_node( ii )->complete_state_assignment();
	}
	update_internal_energy_totals();
	//std::cout << "Set Network State Finished" << std::endl;
	//print_current_state_assignment();
	return total_energy_current_state_assignment_;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyInteractionGraph::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
)
{
	node_considering_alt_state_ = node_ind;

	delta_energy = get_dlazy_node( node_ind )->
		project_deltaE_for_substitution( new_state, prev_energy_for_node );

	total_energy_alternate_state_assignment_ =
		total_energy_current_state_assignment_ + delta_energy;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyInteractionGraph::commit_considered_substitution()
{
	get_dlazy_node( node_considering_alt_state_ )->commit_considered_substitution();

	total_energy_current_state_assignment_ =
		total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if (num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES) {
		update_internal_energy_totals();
	}

	return total_energy_alternate_state_assignment_;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyInteractionGraph::get_energy_current_state_assignment()
{
	//std::cout << "Num rotamer pair energy calculations performed: " << DoubleLazyNode::num_rpe_calcs << std::endl;
	//std::cout << "Num procrastinated comps committed: " << DoubleLazyNode::num_procrastinated_committed << std::endl;
	update_internal_energy_totals();
	return total_energy_current_state_assignment_;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
int
DoubleLazyInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for (std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter) {
		sum += ((DoubleLazyEdge*) *iter)->get_two_body_table_size();
	}
	return sum;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyInteractionGraph::print_current_state_assignment() const
{
	std::cout << "State Assignment: " << std::endl;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		std::cout << "Node " << ii << " state " << get_dlazy_node(ii)->get_current_state() << std::endl;
		get_dlazy_node(ii)->print();
	}

	for (std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter) {
		((DoubleLazyEdge*) (*iter))->print_current_energy();
	}
	std::cout << "Energy: " << total_energy_current_state_assignment_ << std::endl;
}

// @ brief O(1) total energy report.  Protected read access for derived classes.
core::PackerEnergy DoubleLazyInteractionGraph::get_energy_PD_current_state_assignment()
{
	return total_energy_current_state_assignment_;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyInteractionGraph::set_errorfull_deltaE_threshold( core::PackerEnergy )
{}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
core::PackerEnergy
DoubleLazyInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	core::PackerEnergy esum = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_dlazy_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id )) {
			esum += ((DoubleLazyEdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}

	return esum;
}

/// @details DoubleLazyInteractionGraph will return aa submatrices as requested.
bool
DoubleLazyInteractionGraph::aa_submatrix_energies_retrievable() const
{
	return true;
}

int DoubleLazyInteractionGraph::aatype_for_node_state(
	int node_ind,
	int node_state
) const
{
	return get_dlazy_node(node_ind )->aatype_for_state( node_state );
}

ObjexxFCL::FArray2D< core::PackerEnergy >
DoubleLazyInteractionGraph::get_aa_submatrix_energies_for_edge(
	int node1,
	int node2,
	int node1aa,
	int node2aa
) const
{
	return get_dlazy_edge( node1, node2 )->get_aa_submatrix_energies( node1aa, node2aa );
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////

unsigned int
DoubleLazyInteractionGraph::count_static_memory() const
{
	return sizeof( DoubleLazyInteractionGraph );
}

unsigned int
DoubleLazyInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyInteractionGraph::count_dynamic_memory();
	return total_memory;
}

void DoubleLazyInteractionGraph::set_memory_max_for_rpes( int setting )
{
	memory_max_for_rpes_ = setting;
}

void DoubleLazyInteractionGraph::note_submatrix_added(
	int edge_index,
	int submatrix_index, // count from 0 to make module arithmetic easier
	int submatrix_size // in bytes
) const
{
	if ( memory_max_for_rpes_ == 0 ) return;
	curr_memory_for_rpes_ += submatrix_size;
	note_submatrix_accessed( edge_index, submatrix_index );

	while ( curr_memory_for_rpes_ > memory_max_for_rpes_ ) {
		int global_submatrix_index = aa_submatrix_history_list_->tail();
	debug_assert( global_submatrix_index ); /// should never be zero if curr_memroy_for_rpes > memory_max_for_rpes_
		aa_submatrix_history_list_->remove( global_submatrix_index );

		--global_submatrix_index; // convert to a 0-based index for modular arithmetic

		int const edge_index = global_submatrix_index / sqr_num_aa_types_;
		int const local_submatrix_index = global_submatrix_index % sqr_num_aa_types_;

		int savings = dlazy_edge_vector_[ edge_index ]->drop_aa_submatrix( local_submatrix_index );
		//std::cout << "submatrix " << global_submatrix_index << " being deleted " << std::endl;
		curr_memory_for_rpes_ -= savings;
	}
}

void DoubleLazyInteractionGraph::note_submatrix_accessed(
	int edge_index,
	int submatrix_index
) const
{
	if ( memory_max_for_rpes_ == 0 ) return;

	int global_submatrix_index = edge_index * sqr_num_aa_types_ + submatrix_index + 1;
	aa_submatrix_history_list_->move_to_front( global_submatrix_index );
	//std::cout << "submatrix " << global_submatrix_index << " accessed " << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
NodeBase*
DoubleLazyInteractionGraph::create_new_node( int node_index, int num_states)
{
	return new DoubleLazyNode( this, node_index, num_states );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
EdgeBase*
DoubleLazyInteractionGraph::create_new_edge( int index1, int index2)
{
	return new DoubleLazyEdge( this, index1, index2 );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author apl
///
////////////////////////////////////////////////////////////////////////////////
void
DoubleLazyInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		total_energy_current_state_assignment_ += get_dlazy_node( ii )->
			get_one_body_energy_current_state();
	}

	for (std::list<EdgeBase*>::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		total_energy_current_state_assignment_ +=
			((DoubleLazyEdge*) *iter)->get_current_two_body_energy();
	}

	num_commits_since_last_update_ = 0;
	return;
}

} //end namespace pack
}
}
