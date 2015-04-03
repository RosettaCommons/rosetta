// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/LazyInteractionGraph.hh
/// @brief  Interaction graph that computes each rotamer pair energy at most once
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_LazyInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_LazyInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/LazyInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1A.hh>


//Lazy interaction graph uses as much memory as a PDInteractionGraph, but it evaulates
//rotamer pair energies only as it needs them and stores them for later, instead
//of requiring an extensive precomputation phase.

namespace core {
namespace pack {
namespace interaction_graph {

class LazyNode : public OnTheFlyNode
{
public:
	LazyNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states
	);

	virtual ~LazyNode();

	//virtual methods inherited from NodeBase
	virtual void prepare_for_simulated_annealing();
	virtual void print() const;
	virtual bool state_unassigned() const { return current_state_ == 0;}
	virtual core::PackerEnergy get_totalE() const { return curr_state_total_energy_;}

	void assign_zero_state();
	void assign_state(int new_state);
	void partial_assign_state( int new_state );
	void complete_state_assignment();

	inline
	int get_current_state() const { return current_state_; }

	inline
	core::PackerEnergy get_one_body_energy_current_state() const
	{ return curr_state_one_body_energy_; }

	inline
	core::PackerEnergy project_deltaE_for_substitution
	(
		int alternate_state,
		core::PackerEnergy & prev_node_energy
	);

	void commit_considered_substitution();
	core::PackerEnergy compute_pair_energy_for_current_state( int edge_making_energy_request );

	inline
	void acknowledge_neighbors_state_substitution(
		int edge_to_altered_neighbor,
		core::PackerEnergy new_edge_energy,
		int other_node_new_state,
		SparseMatrixIndex const & other_node_new_state_sparse_info
	);

	void
	acknowledge_neighbors_partial_state_substitution(
		int edge_to_altered_neighbor,
		int other_node_new_,
		SparseMatrixIndex const & other_node_new_state_sparse_info
	);

	inline
	SparseMatrixIndex const &
	get_sparse_mat_info_for_curr_state() const;

	void print_internal_energies() const;

	void update_internal_energy_sums();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	virtual
	int aatype_for_state( int state ) const;

protected:

	/*inline
	rotamer_trie const &
	get_current_rotamer()
	{
		return get_rotamer( current_state_ );
	}*/

	/*
	inline
	core::PackerEnergy &
	get_curr_rotamer_actcoords();
	*/

	//Hooks for SASANode< V, E, G > class
	core::PackerEnergy get_curr_pd_energy_total() const { return curr_state_total_energy_;}
	core::PackerEnergy get_alt_pd_energy_total() const { return alternate_state_total_energy_;}
	void set_alternate_state( int alt ) { alternate_state_ = alt; }
	int get_alternate_state() const { return alternate_state_; }
	void calc_deltaEpd( int alternate_state );
	bool considering_alternate_state() const
	{
		return alternate_state_is_being_considered_;
	}


private:

	//void
	//set_alt_aa_offsets_from_edge(int edge_index, FArray2D_int const & offsets);

	void update_internal_vectors();

	/// Pointer Downcasts
	inline LazyEdge const *             get_incident_lazy_edge( int index ) const;
	inline LazyEdge *                   get_incident_lazy_edge( int index )      ;
	inline LazyNode const *             get_adjacent_lazy_node( int index ) const;
	inline LazyNode *                   get_adjacent_lazy_node( int index )      ;
	inline LazyInteractionGraph const * get_lazy_ig_owner() const;
	inline LazyInteractionGraph *       get_lazy_ig_owner()      ;

private:
	// DATA

	ObjexxFCL::FArray3D_int aa_offsets_for_edges_;
	ObjexxFCL::FArray2D_int num_states_for_aa_type_for_higher_indexed_neighbor_;
	std::vector< ObjexxFCL::FArray1A< core::PackerEnergy > > edge_matrix_ptrs_;
	std::vector< int > neighbors_curr_state_;
	std::vector< SparseMatrixIndex > neighbors_curr_state_sparse_info_;

	int current_state_;
	SparseMatrixIndex curr_state_sparse_mat_info_;
	core::PackerEnergy curr_state_one_body_energy_;
	core::PackerEnergy curr_state_total_energy_;
	std::vector< core::PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
	SparseMatrixIndex alt_state_sparse_mat_info_;
	core::PackerEnergy alternate_state_one_body_energy_;
	core::PackerEnergy alternate_state_total_energy_;
	std::vector< core::PackerEnergy > alternate_state_two_body_energies_;

	bool alternate_state_is_being_considered_;
	bool procrastinated_;

	//no default constructor, uncopyable
	LazyNode();
	LazyNode( LazyNode const & );
	LazyNode & operator = ( LazyNode const & );

};

class LazyEdge : public OnTheFlyEdge
{
public:
	LazyEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind
	);

	virtual ~LazyEdge();

	virtual
	void
	set_sparse_aa_info(
		ObjexxFCL::FArray2_bool const &
	);

	virtual
	void force_aa_neighbors( int node1aa, int node2aa);

	virtual
	void force_all_aa_neighbors();

	virtual
	bool
	get_sparse_aa_info(
		int node1aa,
		int node2aa
	) const;

	virtual core::PackerEnergy get_two_body_energy( int const node1state, int const node2state) const;

	//virtual methods inherited from EdgeBase
	virtual void declare_energies_final();
	virtual void prepare_for_simulated_annealing();

	core::PackerEnergy get_current_two_body_energy() const;

	void acknowledge_state_change(
		int node_ind,
		int new_state,
		SparseMatrixIndex const & new_state_sparse_info,
		core::PackerEnergy & new_energy
	);
	void acknowledge_state_zeroed( int node_ind );

	void acknowledge_partial_state_change(
		int node_ind,
		int new_state,
		SparseMatrixIndex const & new_state_sparse_info
	);
	core::PackerEnergy get_energy_following_partial_state_assignment();

	static
	inline
	core::PackerEnergy get_alternate_state_energy_first_node(
		int first_node_alt_state,
		int second_node_orig_state,
		SparseMatrixIndex const & second_node_orig_state_sparse_info,
		int first_node_state_offset_minus_1,
		int second_node_curr_num_states_per_aatype,
		int aa_neighbor_offset,
		ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table
	);

	static
	inline
	core::PackerEnergy get_alternate_state_energy_second_node(
		int first_node_orig_state,
		int second_node_alt_state,
		SparseMatrixIndex const & first_node_orig_state_sparse_info,
		SparseMatrixIndex const & second_node_alternate_state_sparse_info,
		int second_node_alt_state_num_states_per_aatype,
		int aa_neighbor_offset,
		ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table
	);

	static
	inline
	void store_interaction_energy_first_node(
		//int first_node_alt_state,
		//int second_node_orig_state,
		SparseMatrixIndex const & second_node_orig_state_sparse_info,
		int first_node_state_offset_minus_1,
		int second_node_curr_num_states_per_aatype,
		int aa_neighbor_offset,
		ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table,
		core::PackerEnergy interaction_energy
	);


	static
	inline
	void store_interaction_energy_second_node(
		//int first_node_orig_state,
		//int second_node_alt_state,
		SparseMatrixIndex const & first_node_orig_state_sparse_info,
		SparseMatrixIndex const & second_node_alternate_state_sparse_info,
		int second_node_alt_state_num_states_per_aatype,
		int aa_neighbor_offset,
		ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table,
		core::PackerEnergy interaction_energy
	);

	inline
	void acknowledge_substitution(
		int substituted_node_index,
		core::PackerEnergy const curr_state_energy,
		int nodes_new_state,
		SparseMatrixIndex const & nodes_new_state_sparse_info
	);

	void set_edge_weight( Real weight );

	core::PackerEnergy & get_edge_table_ptr();
	int get_two_body_table_size() const;

	void print_current_energy() const;

	bool build_sc_only_rotamer() { return true;}

	ObjexxFCL::FArray2D_int const & get_offsets_for_aatypes( );
	utility::vector1< int > const & get_second_node_num_states_per_aa();

	ObjexxFCL::FArray2D< core::PackerEnergy >
	get_aa_submatrix_energies(
		int node1aa,
		int node2aa
	) const;

	static core::PackerEnergy const NOT_YET_COMPUTED_ENERGY; //an energy lower than any RPE could ever be

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

protected:

	//Hooks for SASAEdge< V, E, G > class
	void declare_energies_final_no_deletion() { LazyEdge::prepare_for_simulated_annealing(); }
	void prepare_for_simulated_annealing_no_deletion() { LazyEdge::prepare_for_simulated_annealing(); }
	bool pd_edge_table_all_zeros() const { return false;}

private:
	void get_energy_for_state_pair(
		int nodes_states[ 2 ],
		SparseMatrixIndex sparse_matrix_indices[ 2 ]
	);

	inline LazyNode const * get_lazy_node( int index ) const;
	inline LazyNode *       get_lazy_node( int index );
	inline LazyInteractionGraph const * get_lazy_ig_owner() const;
	inline LazyInteractionGraph *       get_lazy_ig_owner();

	void
	wipe_two_body_energies_for_node_state( int node, int state );

	mutable AminoAcidNeighborSparseMatrix< core::PackerEnergy > two_body_energies_;
	core::PackerEnergy curr_state_energy_;
	bool partial_state_assignment_;
	bool ran_annealing_since_pair_energy_table_cleared_;

	//no default constructor, uncopyable
	LazyEdge();
	LazyEdge( LazyEdge const & );
	LazyEdge & operator = ( LazyEdge const & );

};

class LazyInteractionGraph : public OnTheFlyInteractionGraph
{
public:
	LazyInteractionGraph( int numNodes );
	virtual ~LazyInteractionGraph();

	//virtual methods inherited from InteractionGraphBase
	virtual void  blanket_assign_state_0();
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state);
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states);
	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	);
	virtual core::PackerEnergy  commit_considered_substitution();
	virtual core::PackerEnergy get_energy_current_state_assignment();
	virtual int get_edge_memory_usage() const;
	virtual void print_current_state_assignment() const;
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE );
	virtual core::PackerEnergy get_energy_sum_for_vertex_group( int group_id );

	/// @brief Override the FixedBBInteractionGraph class's implementation of this function
	/// to return 'true'.
	virtual
	bool
	aa_submatrix_energies_retrievable() const;

	virtual
	int aatype_for_node_state(
		int node_ind,
		int node_state
	) const;

	virtual
	ObjexxFCL::FArray2D< core::PackerEnergy >
	get_aa_submatrix_energies_for_edge(
		int node1,
		int node2,
		int node1aa,
		int node2aa
	) const;

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

protected:

	virtual NodeBase* create_new_node( int node_index, int num_states);
	virtual EdgeBase* create_new_edge( int index1, int index2);

	//Hooks for SASAInterationGraph< V, E, G >
	core::PackerEnergy get_energy_PD_current_state_assignment();
	void update_internal_energy_totals();

	inline
	LazyNode const * get_lazy_node(int index) const
	{
		return static_cast< LazyNode const * > (get_node( index ));
	}

	inline
	LazyNode * get_lazy_node(int index)
	{
		return static_cast< LazyNode * > (get_node( index ));
	}

	inline
	LazyEdge const * get_lazy_edge( int node1, int node2 ) const
	{
		return static_cast< LazyEdge const * > (find_edge( node1, node2 ));
	}

	inline
	LazyEdge * get_lazy_edge( int node1, int node2 )
	{
		return static_cast< LazyEdge * > (find_edge( node1, node2 ));
	}

private:
	int num_aa_types_;
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	LazyInteractionGraph();
	LazyInteractionGraph( LazyInteractionGraph const & );
	LazyInteractionGraph & operator = ( LazyInteractionGraph const & );

};

inline
LazyEdge const * LazyNode::get_incident_lazy_edge( int index ) const
{
	return static_cast< LazyEdge const * >  (get_incident_edge( index ));
}

inline
LazyEdge * LazyNode::get_incident_lazy_edge( int index )
{
	return static_cast< LazyEdge * >  (get_incident_edge( index ));
}

inline
LazyNode const * LazyNode::get_adjacent_lazy_node( int index ) const
{
	return static_cast< LazyNode const * > (get_adjacent_node( index ));
}

inline
LazyNode * LazyNode::get_adjacent_lazy_node( int index )
{
	return static_cast< LazyNode * > (get_adjacent_node( index ));
}

inline
LazyInteractionGraph const * LazyNode::get_lazy_ig_owner() const
{
	return static_cast< LazyInteractionGraph const * > (get_owner());
}

inline
LazyInteractionGraph * LazyNode::get_lazy_ig_owner()
{
	return static_cast< LazyInteractionGraph * > (get_owner());
}


inline
LazyNode const * LazyEdge::get_lazy_node( int index ) const
{
	return static_cast< LazyNode const * > (get_node( index ));
}

inline
LazyNode * LazyEdge::get_lazy_node( int index )
{
	return static_cast< LazyNode * > (get_node( index ));
}

inline
LazyInteractionGraph const * LazyEdge::get_lazy_ig_owner() const
{
	return static_cast< LazyInteractionGraph const * > (get_owner());
}

inline
LazyInteractionGraph * LazyEdge::get_lazy_ig_owner()
{
	return static_cast< LazyInteractionGraph * > (get_owner());
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
inline
void
LazyEdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state,
	SparseMatrixIndex const & nodes_new_state_sparse_info
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = curr_state_energy;

	get_lazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state,
		nodes_new_state_sparse_info
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
inline
void
LazyNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	core::PackerEnergy new_edge_energy,
	int other_node_new_state,
	SparseMatrixIndex const & other_node_new_state_sparse_info
)
{
	curr_state_total_energy_ +=
		new_edge_energy - curr_state_two_body_energies_[edge_to_altered_neighbor];
	curr_state_two_body_energies_[edge_to_altered_neighbor] = new_edge_energy;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	neighbors_curr_state_sparse_info_[ edge_to_altered_neighbor ]  =
		other_node_new_state_sparse_info;
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
inline
void
LazyEdge::store_interaction_energy_first_node(
	SparseMatrixIndex const & second_node_orig_state_sparse_info,
	int first_node_state_offset_minus_1,
	int second_node_curr_num_states_per_aatype,
	int aa_neighbor_offset,
	ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table,
	core::PackerEnergy interaction_energy
)
{

//debug_assert(first_node_alt_state != 0 && second_node_orig_state != 0);
	AminoAcidNeighborSparseMatrix< core::PackerEnergy >::set(
		second_node_orig_state_sparse_info,
		first_node_state_offset_minus_1,
		second_node_curr_num_states_per_aatype,
		aa_neighbor_offset,
		edge_energy_table,
		interaction_energy);
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
inline
void
LazyEdge::store_interaction_energy_second_node(
	SparseMatrixIndex const & first_node_orig_state_sparse_info,
	SparseMatrixIndex const & second_node_alternate_state_sparse_info,
	int second_node_alt_state_num_states_per_aatype,
	int aa_neighbor_offset,
	ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table,
	core::PackerEnergy interaction_energy
)
{
//debug_assert(first_node_orig_state != 0 && second_node_alt_state != 0);
	AminoAcidNeighborSparseMatrix< core::PackerEnergy >::set(
		first_node_orig_state_sparse_info,
		second_node_alternate_state_sparse_info,
		second_node_alt_state_num_states_per_aatype,
		aa_neighbor_offset,
		edge_energy_table,
		interaction_energy);
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
inline
core::PackerEnergy
LazyEdge::get_alternate_state_energy_second_node(
	int first_node_orig_state,
	int second_node_alt_state,
	SparseMatrixIndex const & first_node_orig_state_sparse_info,
	SparseMatrixIndex const & second_node_alternate_state_sparse_info,
	int second_node_alt_state_num_states_per_aatype,
	int aa_neighbor_offset,
	ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table
)
{

	if (first_node_orig_state == 0 || second_node_alt_state == 0) {
		return core::PackerEnergy( 0.0 );
	} else {
		return AminoAcidNeighborSparseMatrix< core::PackerEnergy >::get(
			first_node_orig_state_sparse_info,
			second_node_alternate_state_sparse_info,
			second_node_alt_state_num_states_per_aatype,
			aa_neighbor_offset,
			edge_energy_table );
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
inline
core::PackerEnergy
LazyNode::project_deltaE_for_substitution
(
	int alternate_state,
	core::PackerEnergy & prev_node_energy
)
{
	alternate_state_is_being_considered_ = true;
	//procrastinated_ = false;
	//std::cout << "proj_deltaE: node -  " << get_node_index()
	// << " alt state " << alternate_state << "...";

	alternate_state_ = alternate_state;

	alt_state_sparse_mat_info_ = get_sparse_mat_info_for_state( alternate_state );
	alternate_state_one_body_energy_ = get_one_body_energy( alternate_state );
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
	prev_node_energy = curr_state_total_energy_;

	int aa_neighb_linear_index_offset = aa_offsets_for_edges_.
		index(1, 1, alt_state_sparse_mat_info_.get_aa_type() ) - 1;

	int alt_state_num_states_per_aa_type =
		get_num_states_for_aa_type( alt_state_sparse_mat_info_.get_aa_type() );
	int alt_state_for_aa_type_minus_1 =
		alt_state_sparse_mat_info_.get_state_ind_for_this_aa_type() - 1;
	int nstates_offset =
		num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1) - 1;

	for (int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes();
			++ii, aa_neighb_linear_index_offset += get_num_aa_types()) {

		alternate_state_two_body_energies_[ ii ] =
			get_incident_lazy_edge(ii)->
			get_alternate_state_energy_second_node(
			neighbors_curr_state_[ii],
			alternate_state_,
			neighbors_curr_state_sparse_info_[ ii ],
			alt_state_sparse_mat_info_,
			alt_state_num_states_per_aa_type,
			aa_offsets_for_edges_[
				aa_neighb_linear_index_offset +
				neighbors_curr_state_sparse_info_[ii].get_aa_type()
			],
			edge_matrix_ptrs_[ii]
		);
	}

	for (int ii = get_num_edges_to_smaller_indexed_nodes() + 1;
			ii <= get_num_incident_edges();
			++ii, aa_neighb_linear_index_offset += get_num_aa_types(),
			nstates_offset += get_num_aa_types()) {
		alternate_state_two_body_energies_[ ii ] =
			get_incident_lazy_edge(ii)->
			get_alternate_state_energy_first_node(
			alternate_state_,
			neighbors_curr_state_[ii],
			neighbors_curr_state_sparse_info_[ii],
			alt_state_for_aa_type_minus_1,
			num_states_for_aa_type_for_higher_indexed_neighbor_[
				nstates_offset +
				neighbors_curr_state_sparse_info_[ii].get_aa_type()
			],
			aa_offsets_for_edges_[
				aa_neighb_linear_index_offset +
				neighbors_curr_state_sparse_info_[ii].get_aa_type()
			],
			edge_matrix_ptrs_[ii]
		);
	}

	bool all_energies_computed = true;
	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		if ( alternate_state_two_body_energies_[ ii ] == LazyEdge::NOT_YET_COMPUTED_ENERGY ) {
			all_energies_computed = false;
		}
	}
	if ( all_energies_computed ) {
		return alternate_state_total_energy_ - curr_state_total_energy_;
	}

	aa_neighb_linear_index_offset = aa_offsets_for_edges_.
		index(1, 1, alt_state_sparse_mat_info_.get_aa_type() ) - 1;
	nstates_offset =
		num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1) - 1;

	alternate_state_total_energy_ = alternate_state_one_body_energy_;

	for (int ii = 1; ii <= get_num_incident_edges(); ++ii, aa_neighb_linear_index_offset += get_num_aa_types()) {
		if ( alternate_state_two_body_energies_[ ii ] == LazyEdge::NOT_YET_COMPUTED_ENERGY ) {

			alternate_state_two_body_energies_[ ii ] = compute_rotamer_pair_energy(
				ii,
				alternate_state_,
				neighbors_curr_state_[ ii ]
			);

			if ( ii <= get_num_edges_to_smaller_indexed_nodes() ) {
				LazyEdge::store_interaction_energy_second_node(
					neighbors_curr_state_sparse_info_[ ii ],
					alt_state_sparse_mat_info_,
					alt_state_num_states_per_aa_type,
					aa_offsets_for_edges_[
						aa_neighb_linear_index_offset +
						neighbors_curr_state_sparse_info_[ii].get_aa_type()
					],
					edge_matrix_ptrs_[ii],
					alternate_state_two_body_energies_[ii]
				);
			} else {
				LazyEdge::store_interaction_energy_first_node(
					neighbors_curr_state_sparse_info_[ii],
					alt_state_for_aa_type_minus_1,
					num_states_for_aa_type_for_higher_indexed_neighbor_[
						nstates_offset +
						neighbors_curr_state_sparse_info_[ii].get_aa_type()
					],
					aa_offsets_for_edges_[
						aa_neighb_linear_index_offset +
						neighbors_curr_state_sparse_info_[ii].get_aa_type()
					],
					edge_matrix_ptrs_[ii],
					alternate_state_two_body_energies_[ ii ]
				);
			}
		}
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		if ( ii > get_num_edges_to_smaller_indexed_nodes() ) {
			nstates_offset += get_num_aa_types();
		}
	}

	//std::cout << " " << (double) num_rpe_calcs_this_sub / get_num_incident_edges();

	return alternate_state_total_energy_ - curr_state_total_energy_;

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
inline
core::PackerEnergy
LazyEdge::get_alternate_state_energy_first_node(
	int first_node_alt_state,
	int second_node_orig_state,
	SparseMatrixIndex const & second_node_orig_state_sparse_info,
	int first_node_state_offset_minus_1,
	int second_node_curr_num_states_per_aatype,
	int aa_neighbor_offset,
	ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table

)
{
	if ( first_node_alt_state == 0 || second_node_orig_state == 0 ) {
		return core::PackerEnergy( 0.0 );
	} else {
		return AminoAcidNeighborSparseMatrix< core::PackerEnergy >::get(
			second_node_orig_state_sparse_info,
			first_node_state_offset_minus_1,
			second_node_curr_num_states_per_aatype,
			aa_neighbor_offset,
			edge_energy_table );
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
inline
SparseMatrixIndex const &
LazyNode::get_sparse_mat_info_for_curr_state() const
{
	return get_sparse_mat_info_for_state( current_state_ );
}

}
}
}

#endif
