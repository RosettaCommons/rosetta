// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/LazyInteractionGraph.cc
/// @brief  Interaction graph that computes each rotamer pair energy at most once
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <core/pack/interaction_graph/LazyInteractionGraph.hh>

/// C++ headers
#include <iostream>

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
LazyNode::LazyNode(
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
LazyNode::~LazyNode()
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
void
LazyNode::prepare_for_simulated_annealing()
{
	if ( ! get_edge_vector_up_to_date() ) update_internal_vectors();
	/*for (int ii = 1; ii <= get_num_states(); ++ii) {
	mark_coordinates_current( ii ); /// What did this used to to?
	}*/
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
LazyNode::print() const
{

	std::cout << "LazyNode " << get_node_index() << " with " << get_num_states() << " states" << std::endl;
	std::cout << "curr_state " << current_state_ << " ";
	std::cout << "curr_state_sparse_mat_info_ ";
	std::cout << curr_state_sparse_mat_info_.get_aa_type() << " ";
	std::cout << curr_state_sparse_mat_info_.get_state_ind_for_this_aa_type() << " ";
	std::cout << "Curr One Body Energy: " << curr_state_one_body_energy_ << std::endl;
	std::cout << "Curr Two Body Energies:";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		std::cout << " " << get_index_of_adjacent_node(ii) << ":" << curr_state_two_body_energies_[ ii ];
	}
	std::cout << std::endl;

	if ( ! alternate_state_is_being_considered_ ) return;
	std::cout << "Alt One Body Energy: " << alternate_state_one_body_energy_ << std::endl;
	std::cout << "Alt Two Body Energies:";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
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
LazyNode::count_static_memory() const
{
	return sizeof( LazyNode );
}

unsigned int
LazyNode::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyNode::count_dynamic_memory();

	total_memory += neighbors_curr_state_.size() * sizeof (int );

	total_memory += aa_offsets_for_edges_.size() * sizeof( int );
	total_memory += num_states_for_aa_type_for_higher_indexed_neighbor_.size() * sizeof( int );
	total_memory += neighbors_curr_state_.size() * sizeof( int );
	total_memory += neighbors_curr_state_sparse_info_.size() * sizeof( SparseMatrixIndex );
	total_memory += edge_matrix_ptrs_.size() * sizeof( ObjexxFCL::FArray1A< core::PackerEnergy > );

	total_memory += curr_state_two_body_energies_.size() * sizeof( core::PackerEnergy );
	total_memory += alternate_state_two_body_energies_.size() * sizeof( core::PackerEnergy );

	return total_memory;
}

int
LazyNode::aatype_for_state( int state ) const
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
LazyNode::assign_zero_state()
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

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_lazy_edge(ii)->acknowledge_state_zeroed( get_node_index() );
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
LazyNode::assign_state(int new_state)
{
	debug_assert( new_state >= 0 && new_state <= get_num_states());

	if ( new_state == 0 ) {
		assign_zero_state();
	} else {
		//std::cout << "assign_state: node -  " << get_node_index() <<
		// " new state " << new_state << "...";
		current_state_ = new_state;
		curr_state_sparse_mat_info_ = get_sparse_mat_info_for_state( current_state_ );
		curr_state_one_body_energy_ = get_one_body_energy( current_state_ );
		curr_state_total_energy_ = curr_state_one_body_energy_;
		alternate_state_is_being_considered_ = false;

		for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
			get_incident_lazy_edge(ii)->acknowledge_state_change(
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
LazyNode::partial_assign_state( int new_state )
{
	if ( new_state == 0 ) {
		assign_zero_state();
		return;
	}

	current_state_ = new_state;
	curr_state_sparse_mat_info_ =
		get_sparse_mat_info_for_state( current_state_ );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_lazy_edge(ii)->acknowledge_partial_state_change(
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
void LazyNode::complete_state_assignment()
{
	if ( current_state_ == 0 ) return;

	curr_state_total_energy_ = curr_state_one_body_energy_ =
		get_one_body_energy( current_state_ );
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_two_body_energies_[ ii ] =
			get_incident_lazy_edge( ii )->
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
LazyNode::commit_considered_substitution()
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
	// curr_state_total_energy_ = curr_state_one_body_energy_;
	// for (int ii = 1; ii <= get_num_incident_edges(); ++ii)
	// {
	//  if ( curr_state_two_body_energies_[ ii ] == LazyEdge::NOT_YET_COMPUTED_ENERGY )
	//  {
	//   curr_state_two_body_energies_[ ii ] = compute_pair_energy_for_current_state( ii );
	//  }
	//  curr_state_total_energy_ += curr_state_two_body_energies_[ ii ];
	// }
	// procrastinated_ = false;
	// ++num_procrastinated_committed;
	//}


	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		get_incident_lazy_edge(ii)->acknowledge_substitution(
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
LazyNode::compute_pair_energy_for_current_state(
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
LazyNode::acknowledge_neighbors_partial_state_substitution(
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
LazyNode::print_internal_energies() const
{
	std::cout << "curr_state " << current_state_ << " ";
	std::cout << "curr_state_sparse_mat_info_ ";
	std::cout << curr_state_sparse_mat_info_.get_aa_type() << " ";
	std::cout << curr_state_sparse_mat_info_.get_state_ind_for_this_aa_type() << " ";
	std::cout << "curr_state_one_body_energy_ ";
	std::cout << curr_state_one_body_energy_ << " ";
	std::cout << "curr_state_total_energy_" << curr_state_total_energy_ << " ";
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
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
LazyNode::update_internal_energy_sums()
{
	debug_assert( get_edge_vector_up_to_date() );
	curr_state_total_energy_ = 0;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		curr_state_total_energy_ += get_incident_lazy_edge(ii)->get_current_two_body_energy();
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
void LazyNode::update_internal_vectors()
{
	NodeBase::update_edge_vector();
	neighbors_curr_state_.resize( get_num_incident_edges() + 1);
	neighbors_curr_state_sparse_info_.resize( get_num_incident_edges() + 1);

	edge_matrix_ptrs_.clear();
	edge_matrix_ptrs_.reserve( get_num_incident_edges() + 1);
	edge_matrix_ptrs_.push_back( ObjexxFCL::FArray1A< core::PackerEnergy >() ); //occupy the 0th position

	aa_offsets_for_edges_.dimension(
		get_num_aa_types(), get_num_incident_edges(), get_num_aa_types());
	num_states_for_aa_type_for_higher_indexed_neighbor_.dimension(
		get_num_aa_types(), get_num_edges_to_larger_indexed_nodes());

	//copy offsets from edges
	//int neighb_aa_offset =
	// num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1);
	int count_neighbs_with_higher_indices = 0;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		neighbors_curr_state_sparse_info_[ii].set_aa_type( 1 );

		core::PackerEnergy & edge_table_ref =
			get_incident_lazy_edge(ii)->get_edge_table_ptr();
		edge_matrix_ptrs_.push_back( ObjexxFCL::FArray1A< core::PackerEnergy >( edge_table_ref ));


		int edge_table_size = get_incident_lazy_edge(ii)->get_two_body_table_size();
		edge_matrix_ptrs_[ii].dimension( edge_table_size );

		ObjexxFCL::FArray2D_int const & edge_aa_neighb_offsets =
			get_incident_lazy_edge(ii)->get_offsets_for_aatypes();
		utility::vector1< int > const & neighb_num_states_per_aa =
			get_incident_lazy_edge(ii)->get_second_node_num_states_per_aa();

		if ( get_node_index() < get_index_of_adjacent_node(ii) ) {
			++count_neighbs_with_higher_indices;
			for ( int jj = 1; jj <= get_num_aa_types(); ++jj ) {
				for ( int kk = 1; kk <= get_num_aa_types(); ++kk ) {
					aa_offsets_for_edges_(kk, ii, jj) = edge_aa_neighb_offsets(kk, jj);
				}
				num_states_for_aa_type_for_higher_indexed_neighbor_(
					jj, count_neighbs_with_higher_indices) =
					neighb_num_states_per_aa[ jj ];
				//++neighb_aa_offset;
			}
		} else {
			for ( int jj = 1; jj <= get_num_aa_types(); ++jj ) {
				for ( int kk = 1; kk <= get_num_aa_types(); ++kk ) {
					aa_offsets_for_edges_(kk, ii, jj) = edge_aa_neighb_offsets(jj, kk);
				}
			}
		}
	}

	curr_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	alternate_state_two_body_energies_.resize( get_num_incident_edges() + 1);
	return;
}

// @ brief - allow derived class to "drive" through the deltaE calculation
void
LazyNode::calc_deltaEpd( int alternate_state )
{
	core::PackerEnergy dummy(0.0f);
	project_deltaE_for_substitution( alternate_state, dummy );
}


//-----------------------------------------------------------------//

core::PackerEnergy const LazyEdge::NOT_YET_COMPUTED_ENERGY = -1234;

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
LazyEdge::LazyEdge(
	InteractionGraphBase* owner,
	int first_node_ind,
	int second_node_ind
):
	OnTheFlyEdge( owner, first_node_ind, second_node_ind),
	//energy_table_size_(0)
	two_body_energies_(
	get_lazy_node(0)->get_num_states_for_aa_types(),
	get_lazy_node(1)->get_num_states_for_aa_types()
	),
	curr_state_energy_( 0.0f ),
	partial_state_assignment_( false ),
	ran_annealing_since_pair_energy_table_cleared_( false )
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
LazyEdge::~LazyEdge()
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
void
LazyEdge::set_sparse_aa_info(
	ObjexxFCL::FArray2_bool const & aa_neighbors
)
{
	two_body_energies_.set_sparse_aa_info( aa_neighbors );
	two_body_energies_.blanket_set( NOT_YET_COMPUTED_ENERGY );
	ran_annealing_since_pair_energy_table_cleared_ = false;
}

/// @brief returns whether two amino acid types are represented as neighbors
bool LazyEdge::get_sparse_aa_info( int node1aa, int node2aa ) const
{
	return two_body_energies_.get_sparse_aa_info( node1aa, node2aa );
}

/// @brief re-allocates two-body energy table after forcing a pair of amino acids
/// to become neighbors that were not initially declared to be neighbors
///
/// @param node1aa - [in] - the amino acid type for the node with the smaller index
/// @param node2aa - [in] - the amino acid type for the node with the larger index
///
void LazyEdge::force_aa_neighbors(int node1aa, int node2aa)
{
	two_body_energies_.force_aa_neighbors( node1aa, node2aa );
	two_body_energies_.blanket_set( NOT_YET_COMPUTED_ENERGY );
	ran_annealing_since_pair_energy_table_cleared_ = false;
}

/// @brief re-allocates two-body energy table after forcing a pair of amino acids
/// to become neighbors that were not initially declared to be neighbors
///
/// @param node1aa - [in] - the amino acid type for the node with the smaller index
/// @param node2aa - [in] - the amino acid type for the node with the larger index
///
void LazyEdge::force_all_aa_neighbors()
{
	two_body_energies_.force_all_aa_neighbors();
	two_body_energies_.blanket_set( NOT_YET_COMPUTED_ENERGY );
	ran_annealing_since_pair_energy_table_cleared_ = false;
}


core::PackerEnergy LazyEdge::get_two_body_energy(
	int const node1state,
	int const node2state
) const
{
	SparseMatrixIndex node1info = get_otf_node(0)->get_sparse_mat_info_for_state( node1state );
	SparseMatrixIndex node2info = get_otf_node(1)->get_sparse_mat_info_for_state( node2state );

	core::PackerEnergy energy = two_body_energies_.get( node1info, node2info );

	if ( energy == NOT_YET_COMPUTED_ENERGY ) {
		energy = get_otf_node(0)->compute_rotamer_pair_energy(
			get_edges_position_in_nodes_edge_vector( 0 ),
			node1state, node2state );
		two_body_energies_.set( node1info, node2info, energy );
	}
	return energy;
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
LazyEdge::declare_energies_final()
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
void
LazyEdge::prepare_for_simulated_annealing()
{
	if ( ! ran_annealing_since_pair_energy_table_cleared_ ) {
		ran_annealing_since_pair_energy_table_cleared_ = true;
		if ( two_body_energies_.size() == 0 ) {
			delete this;
			return;
		}
		return;
	} else {
		/*for (int ii = 0; ii < 2; ++ii) {
		if ( get_lazy_node( ii )->get_any_coordinates_not_current() ) {
		for (int jj = 1; jj <= get_lazy_node( ii )->get_num_states(); ++jj ) {
		if ( ! get_lazy_node( ii )->get_coordinates_current( jj ) ) {
		wipe_two_body_energies_for_node_state( ii, jj );
		}
		}
		}
		}*/
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

unsigned int
LazyEdge::count_static_memory() const
{
	return sizeof( LazyEdge );
}


unsigned int
LazyEdge::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyEdge::count_dynamic_memory();
	total_memory += two_body_energies_.get_table_size() * sizeof( core::PackerEnergy );
	total_memory += two_body_energies_.get_offset_table_size_in_bytes();

	return total_memory;
}


void
LazyEdge::set_edge_weight( Real weight )
{
	Real const scale_factor = weight / edge_weight();
	for ( int ii = 1; ii <= two_body_energies_.size(); ++ii ) {
		if ( two_body_energies_[ ii ] != NOT_YET_COMPUTED_ENERGY ) {
			two_body_energies_[ ii ] *= scale_factor;
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
LazyEdge::get_current_two_body_energy() const
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
LazyEdge::acknowledge_state_change(
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
		get_lazy_node( node_not_substituted )->get_current_state();
	nodes_curr_states_sparse_info[ node_not_substituted ] =
		get_lazy_node( node_not_substituted )->
		get_sparse_mat_info_for_curr_state();

	get_energy_for_state_pair( nodes_curr_states, nodes_curr_states_sparse_info );
	new_energy = curr_state_energy_;

	get_lazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
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
LazyEdge::acknowledge_state_zeroed( int node_ind )
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;
	SparseMatrixIndex dummy_sparse_info;
	dummy_sparse_info.set_aa_type( 1 );
	dummy_sparse_info.set_state_ind_for_this_aa_type(1);

	get_lazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
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
void LazyEdge::acknowledge_partial_state_change(
	int node_ind,
	int new_state,
	SparseMatrixIndex const & new_state_sparse_info
)
{
	int node_substituted =  ( node_ind == get_node_index(0) ? 0 : 1);
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = 0;

	get_lazy_node( node_not_substituted )->acknowledge_neighbors_partial_state_substitution(
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
LazyEdge::get_energy_following_partial_state_assignment()
{
	if ( partial_state_assignment_ ) {
		int nodes_curr_states[2];
		SparseMatrixIndex nodes_curr_states_sparse_info[2];

		for ( int ii = 0; ii < 2; ++ii ) {
			nodes_curr_states[ ii ] =
				get_lazy_node( ii )->get_current_state();
			nodes_curr_states_sparse_info[ ii ] =
				get_lazy_node( ii )->get_sparse_mat_info_for_curr_state();
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
core::PackerEnergy & LazyEdge::get_edge_table_ptr()
{
	return two_body_energies_.getMatrixPointer();
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
int LazyEdge::get_two_body_table_size() const
{
	return two_body_energies_.get_table_size();
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
ObjexxFCL::FArray2D_int const &
LazyEdge::get_offsets_for_aatypes( )
{
	return two_body_energies_.getAANeighborOffsets();
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
utility::vector1< int > const &
LazyEdge::get_second_node_num_states_per_aa()
{
	return get_lazy_node(1)->get_num_states_for_aa_types();
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
LazyEdge::get_energy_for_state_pair(
	int nodes_states[ 2 ],
	SparseMatrixIndex sparse_matrix_indices_[ 2 ]
)
{
	bool one_node_in_zero_state = ( nodes_states[0] == 0 || nodes_states[1] == 0 );

	if (  one_node_in_zero_state ) {
		curr_state_energy_ = 0;
	} else {
		curr_state_energy_ = two_body_energies_.get(
			sparse_matrix_indices_[0],
			sparse_matrix_indices_[1]);
		if ( curr_state_energy_ == NOT_YET_COMPUTED_ENERGY ) {
			curr_state_energy_ = get_lazy_node( 0 )->
				compute_pair_energy_for_current_state(
				get_edges_position_in_nodes_edge_vector( 0 ) );
			two_body_energies_.set(
				sparse_matrix_indices_[0],
				sparse_matrix_indices_[1],
				curr_state_energy_);
		}
	}
}

ObjexxFCL::FArray2D< core::PackerEnergy >
LazyEdge::get_aa_submatrix_energies(
	int node1aa,
	int node2aa
) const
{
	ObjexxFCL::FArray2D< core::PackerEnergy > submat = two_body_energies_.get_aa_submatrix_energies( node1aa, node2aa );
	int const iioffset = get_lazy_node(0)->get_state_offset_for_aatype( node1aa );
	int const jjoffset = get_lazy_node(1)->get_state_offset_for_aatype( node2aa );

	for ( Size ii = 1; ii <= submat.size2(); ++ii ) {
		for ( Size jj = 1; jj <= submat.size1(); ++jj ) {
			if ( submat(jj,ii) == NOT_YET_COMPUTED_ENERGY ) {
				int const ii_state = ii + iioffset;
				int const jj_state = jj + jjoffset;
				core::PackerEnergy iijjenergy = get_otf_node(0)->compute_rotamer_pair_energy(
					get_edges_position_in_nodes_edge_vector( 0 ),
					ii_state, jj_state );
				two_body_energies_.set(
					get_otf_node(0)->get_sparse_mat_info_for_state( ii_state ),
					get_otf_node(1)->get_sparse_mat_info_for_state( jj_state ),
					iijjenergy );
				submat(jj,ii) = iijjenergy;
			}
		}
	}
	return submat;
}


void
LazyEdge::wipe_two_body_energies_for_node_state(
	int node,
	int state
)
{
	int other_node = node == 0 ? 1 : 0;

	SparseMatrixIndex states[ 2 ];
	states[ node ] = get_lazy_node( node )
		->get_sparse_mat_info_for_state( state );
	int const other_node_num_states =
		get_lazy_node( other_node )->get_num_states();

	for ( int ii = 1; ii <= other_node_num_states; ++ii ) {
		states[ other_node ] = get_lazy_node( other_node )
			->get_sparse_mat_info_for_state( ii );
		two_body_energies_.set( states[ 0 ], states[ 1 ], NOT_YET_COMPUTED_ENERGY );

	}

}

void
LazyEdge::print_current_energy() const
{
	std::cout << "LazyEdge: " << get_node_index( 0 ) << "/" << get_node_index( 1 );
	std::cout << " energy= " << curr_state_energy_ << std::endl;
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
LazyInteractionGraph::LazyInteractionGraph(
	int numNodes
) : OnTheFlyInteractionGraph( numNodes )
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
LazyInteractionGraph::~LazyInteractionGraph()
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
void
LazyInteractionGraph::blanket_assign_state_0()
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_lazy_node( ii )->assign_zero_state();
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
LazyInteractionGraph::set_state_for_node(int node_ind, int new_state)
{
	get_lazy_node( node_ind )->assign_state( new_state );
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
LazyInteractionGraph::set_network_state( ObjexxFCL::FArray1_int & node_states)
{
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_lazy_node( ii )->partial_assign_state( node_states( ii ) );
	}
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		get_lazy_node( ii )->complete_state_assignment();
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
LazyInteractionGraph::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
)
{
	node_considering_alt_state_ = node_ind;

	delta_energy = get_lazy_node( node_ind )->
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
LazyInteractionGraph::commit_considered_substitution()
{
	get_lazy_node( node_considering_alt_state_ )->commit_considered_substitution();

	total_energy_current_state_assignment_ =
		total_energy_alternate_state_assignment_;

	++num_commits_since_last_update_;
	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
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
LazyInteractionGraph::get_energy_current_state_assignment()
{
	//std::cout << "Num rotamer pair energy calculations performed: " << LazyNode::num_rpe_calcs << std::endl;
	//std::cout << "Num procrastinated comps committed: " << LazyNode::num_procrastinated_committed << std::endl;
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
LazyInteractionGraph::get_edge_memory_usage() const
{
	int sum = 0;
	for ( std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		sum += ((LazyEdge*) *iter)->get_two_body_table_size();
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
LazyInteractionGraph::print_current_state_assignment() const
{
	std::cout << "State Assignment: " << std::endl;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		std::cout << "Node " << ii << " state " << get_lazy_node(ii)->get_current_state() << std::endl;
		get_lazy_node(ii)->print();
	}

	for ( std::list< EdgeBase* >::const_iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		((LazyEdge*) (*iter))->print_current_energy();
	}
	std::cout << "Energy: " << total_energy_current_state_assignment_ << std::endl;
}

// @ brief O(1) total energy report.  Protected read access for derived classes.
core::PackerEnergy LazyInteractionGraph::get_energy_PD_current_state_assignment()
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
LazyInteractionGraph::set_errorfull_deltaE_threshold( core::PackerEnergy )
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
LazyInteractionGraph::get_energy_sum_for_vertex_group( int group_id )
{
	core::PackerEnergy esum = 0;
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		if ( get_vertex_member_of_energy_sum_group( ii, group_id ) ) {
			esum += get_lazy_node( ii )->get_one_body_energy_current_state();
		}
	}

	for ( std::list< EdgeBase* >::iterator edge_iter = get_edge_list_begin();
			edge_iter != get_edge_list_end(); ++edge_iter ) {
		int first_node_ind = (*edge_iter)->get_first_node_ind();
		int second_node_ind = (*edge_iter)->get_second_node_ind();

		if ( get_vertex_member_of_energy_sum_group( first_node_ind, group_id )
				&& get_vertex_member_of_energy_sum_group( second_node_ind, group_id ) ) {
			esum += ((LazyEdge*) (*edge_iter))->get_current_two_body_energy();
		}
	}

	return esum;
}

/// @details LazyInteractionGraph will return aa submatrices as requested.
bool
LazyInteractionGraph::aa_submatrix_energies_retrievable() const
{
	return true;
}

int LazyInteractionGraph::aatype_for_node_state(
	int node_ind,
	int node_state
) const
{
	return get_lazy_node(node_ind )->aatype_for_state( node_state );
}

ObjexxFCL::FArray2D< core::PackerEnergy >
LazyInteractionGraph::get_aa_submatrix_energies_for_edge(
	int node1,
	int node2,
	int node1aa,
	int node2aa
) const
{
	return get_lazy_edge( node1, node2 )->get_aa_submatrix_energies( node1aa, node2aa );
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
LazyInteractionGraph::count_static_memory() const
{
	return sizeof( LazyInteractionGraph );
}

unsigned int
LazyInteractionGraph::count_dynamic_memory() const
{
	unsigned int total_memory = OnTheFlyInteractionGraph::count_dynamic_memory();
	return total_memory;
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
LazyInteractionGraph::create_new_node( int node_index, int num_states)
{
	return new LazyNode( this, node_index, num_states );
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
LazyInteractionGraph::create_new_edge( int index1, int index2)
{
	return new LazyEdge( this, index1, index2 );
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
LazyInteractionGraph::update_internal_energy_totals()
{
	total_energy_current_state_assignment_ = 0;

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		total_energy_current_state_assignment_ += get_lazy_node( ii )->
			get_one_body_energy_current_state();
	}

	for ( std::list<EdgeBase*>::iterator iter = get_edge_list_begin();
			iter != get_edge_list_end(); ++iter ) {
		total_energy_current_state_assignment_ +=
			((LazyEdge*) *iter)->get_current_two_body_energy();
	}

	num_commits_since_last_update_ = 0;
	return;
}

} //end namespace pack
}
}
