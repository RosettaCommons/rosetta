// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/DoubleLazyInteractionGraph.hh
/// @brief  Interaction graph that computes each rotamer pair energy at most once
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_DoubleLazyInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_DoubleLazyInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/DoubleLazyInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>
// AUTO-REMOVED #include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>

#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>

// Utility headers
#include <utility/in_place_list.fwd.hh>
#include <utility/vector0.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.hh>


// The double-lazy interaction graph is lazy both in its computation of
// rotamer pair energies and in its allocation of memory for rotamer
// pair energies.  Each rotamer pair energy is computed at most once.
// AA neighbor submatrices are only allocated
// when they are needed.  This interaction graph is particularly well
// suited for multistate design where single aa-submatrices are used
// at a time.

namespace core {
namespace pack {
namespace interaction_graph {

class DoubleLazyNode : public OnTheFlyNode
{
public:
	DoubleLazyNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states
	);

	virtual ~DoubleLazyNode();

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
	//unsigned int potential_edge_memory_use() const;

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
	inline DoubleLazyEdge const *             get_incident_dlazy_edge( int index ) const;
	inline DoubleLazyEdge *                   get_incident_dlazy_edge( int index )      ;
	inline DoubleLazyNode const *             get_adjacent_dlazy_node( int index ) const;
	inline DoubleLazyNode *                   get_adjacent_dlazy_node( int index )      ;
	inline DoubleLazyInteractionGraph const * get_dlazy_ig_owner() const;
	inline DoubleLazyInteractionGraph *       get_dlazy_ig_owner()      ;

private:
	// DATA

	//FArray3D_int aa_offsets_for_edges_;
	//FArray2D_int num_states_for_aa_type_for_higher_indexed_neighbor_;
	//std::vector< FArray1Da< core::PackerEnergy > > edge_matrix_ptrs_;

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
	DoubleLazyNode();
	DoubleLazyNode( DoubleLazyNode const & );
	DoubleLazyNode & operator = ( DoubleLazyNode const & );

};

class DoubleLazyEdge : public OnTheFlyEdge
{
public:
	DoubleLazyEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind
	);

	virtual ~DoubleLazyEdge();

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

	core::PackerEnergy
	get_two_body_energy_smi(
		int const node1state,
		int const node2state,
		SparseMatrixIndex const & node1info,
		SparseMatrixIndex const & node2info
	) const;

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

	/* inline
	core::PackerEnergy get_alternate_state_energy_first_node(
		int first_node_alt_state,
		int second_node_orig_state,
		SparseMatrixIndex const & first_node_alt_state_sparse_info,
		SparseMatrixIndex const & second_node_orig_state_sparse_info
	); */

	/* inline
	core::PackerEnergy get_alternate_state_energy_second_node(
		int first_node_orig_state,
		int second_node_alt_state,
		SparseMatrixIndex const & first_node_orig_state_sparse_info,
		SparseMatrixIndex const & second_node_alternate_state_sparse_info
	); */

	inline
	void acknowledge_substitution(
		int substituted_node_index,
		core::PackerEnergy const curr_state_energy,
		int nodes_new_state,
		SparseMatrixIndex const & nodes_new_state_sparse_info
	);

	void set_edge_weight( Real weight );

	//core::PackerEnergy & get_edge_table_ptr();
	int get_two_body_table_size() const;

	void print_current_energy() const;

	bool build_sc_only_rotamer() { return true;}

	//FArray2D_int const & get_offsets_for_aatypes( );
	//utility::vector1< int > const & get_second_node_num_states_per_aa();

	ObjexxFCL::FArray2D< core::PackerEnergy >
	get_aa_submatrix_energies(
		int node1aa,
		int node2aa
	) const;

	static core::PackerEnergy const NOT_YET_COMPUTED_ENERGY; //an energy lower than any RPE could ever be

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	unsigned int
	actual_twobody_memory_use() const;

	unsigned int
	potential_twobody_memory_use() const;

	/// @brief For use by the DoubleLazyInteractionGraph only;
	/// sets the edge index for this edge -- this is an arbitrary
	/// index, but is used to so that the DLIG can communicate back
	/// and forth with this edge.
	void
	set_edge_index(
		int index
	);

	/// @brief For use by the DoubleLazyInteractionGraph only; if the DLIG is in "memory conservation"
	/// mode, then it may request that its edges drop certain amino acid submatrices which have not been activated
	/// for a "long time."  This function will cause an amino-acid submatrix to be dropped, and will return
	/// the number of bytes that have been freed in the process. "const" but accesses mutable data
	int
	drop_aa_submatrix(
		int submat_ind
	) const;

protected:

	//Hooks for SASAEdge< V, E, G > class
	void declare_energies_final_no_deletion() { DoubleLazyEdge::prepare_for_simulated_annealing(); }
	void prepare_for_simulated_annealing_no_deletion() { DoubleLazyEdge::prepare_for_simulated_annealing(); }
	bool pd_edge_table_all_zeros() const { return false;}

private:

	void
	prep_aa_submatrix(
		int node1aa,
		int node2aa
	) const;

	core::PackerEnergy
	read_aa_submatrix(
		SparseMatrixIndex node1info,
		SparseMatrixIndex node2info
	) const;

	void
	set_aa_submatrix(
		SparseMatrixIndex node1info,
		SparseMatrixIndex node2info,
		core::PackerEnergy setting
	) const;

	void get_energy_for_state_pair(
		int nodes_states[ 2 ],
		SparseMatrixIndex sparse_matrix_indices[ 2 ]
	);

	inline DoubleLazyNode const * get_dlazy_node( int index ) const;
	inline DoubleLazyNode *       get_dlazy_node( int index );
	inline DoubleLazyInteractionGraph const * get_dlazy_ig_owner() const;
	inline DoubleLazyInteractionGraph *       get_dlazy_ig_owner();

	void
	wipe_two_body_energies_for_node_state( int node, int state );

	int
	submatrix_index( int nod1aa, int node2aa ) const;

	int
	submatrix_size( int node1aa, int node2aa ) const;

	std::pair< int, int >
	aa_indices_from_submatrix_index( int submat_ind ) const;

private:
	/// record which amino acid pairs are neighbors
	ObjexxFCL::FArray2D< unsigned char > sparse_aa_neighbors_;
	/// Allocate the amino-acid submatrices of the two-body energy table only as needed.
	/// In a more advanced version of this graph, these submatrices would be deleted again
	/// (and possibly re-instantiated!) after a long time had passed since they were last useful.
	mutable ObjexxFCL::FArray2D< ObjexxFCL::FArray2D< core::PackerEnergy > * > two_body_energies_;
	core::PackerEnergy curr_state_energy_;
	bool partial_state_assignment_;
	bool ran_annealing_since_pair_energy_table_cleared_;

	// The double-lazy interaction graph keeps a vector of edges and assigns an
	// index to each edge.  The DLIG and its edges communicate to each other based
	// on this (arbitrary) indexing.  The index is assigned in
	// prepare_for_simulated_annealing
	int edge_index_;

	//no default constructor, uncopyable
	DoubleLazyEdge();
	DoubleLazyEdge( DoubleLazyEdge const & );
	DoubleLazyEdge & operator = ( DoubleLazyEdge const & );

};

/// @brief The double lazy interaction graph is primarily useful for multistate design
/// where one is interested in knowing at particular edge, all of the rotamer pair energies
/// for a particular amino acid assignment.
/// The double lazy interaction graph is lazy in two ways: first, in delaying the computation
/// of rotamer pair energies until they are needed, and second, in delaying the allocation of
/// memory for rotamer pair energies until that memory is needed.
/// The DLIG will do one of two things once it allocates space for a block of rotamer pairs:
/// 1) In its standard operating behavior, it will leave that space allocated until
/// the graph is destroyed, which means that the energies it stores in that block will never
/// be computed more than once; or
/// 2) In its alternate operating behavior, the LMIG will deallocate some of those blocks
/// to make sure that it never uses more than some maximum amount of memory on RPEs.
class DoubleLazyInteractionGraph : public OnTheFlyInteractionGraph
{
public:
	typedef OnTheFlyInteractionGraph parent;
	typedef utility::pointer::owning_ptr< utility::in_place_list< int > > InPlaceIntListOP;

public:
	DoubleLazyInteractionGraph( int numNodes );
	virtual ~DoubleLazyInteractionGraph();
	virtual void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

	//virtual methods inherited from InteractionGraphBase
	virtual void prepare_for_simulated_annealing();
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

	void set_memory_max_for_rpes( int setting );

	/// @brief For use only from the DoubleLazyEdge; the DLE will report
	/// to the DLIG after it adds a submatrix of rotamer pair energies
	/// (that is, all the RPEs for a particular pair of amino acids)
	/// that it has done so; the DLIG monitors the access and the memory
	/// usage of the various submatrices in the graph, and may, during
	/// this function call, request that various edges (possibly the one
	/// invoking this function!) drop some of their submatrices.
	void note_submatrix_added(
		int edge_index,
		int submatrix_index,
		int submatrix_size // in bytes
	) const;

	/// @brief For use only from the DoubleLazyEdge; the DLE will
	/// report to the DLIG after it a rotamer pair energy submatrix
	/// has been read from.  The DLIG keeps track of how recently each
	/// submatrix has been accessed to ensure that, when it does ask
	/// an edge to drop a submatrix, the submatrix being dropped
	/// is the one that was accessed most distantly in the past.
	void note_submatrix_accessed(
		int edge_index,
		int submatrix_index
	) const;


protected:

	virtual NodeBase* create_new_node( int node_index, int num_states);
	virtual EdgeBase* create_new_edge( int index1, int index2);

	//Hooks for SASAInterationGraph< V, E, G >
	core::PackerEnergy get_energy_PD_current_state_assignment();
	void update_internal_energy_totals();

	inline
	DoubleLazyNode const * get_dlazy_node(int index) const
	{
		return static_cast< DoubleLazyNode const * > (get_node( index ));
	}

	inline
	DoubleLazyNode * get_dlazy_node(int index)
	{
		return static_cast< DoubleLazyNode * > (get_node( index ));
	}

	inline
	DoubleLazyEdge const * get_dlazy_edge( int node1, int node2 ) const
	{
		return static_cast< DoubleLazyEdge const * > (find_edge( node1, node2 ));
	}

	inline
	DoubleLazyEdge * get_dlazy_edge( int node1, int node2 )
	{
		return static_cast< DoubleLazyEdge * > (find_edge( node1, node2 ));
	}

private:
	int sqr_num_aa_types_; // == num_aa_types_ * num_aa_types_;
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	// the maximum amount of memory that should be dedicated
	// toward rotamer pair energies within the graph.
	int memory_max_for_rpes_;
	mutable int curr_memory_for_rpes_;

	/// initialized in prepare_for_simulated_annealing();
	utility::vector0< DoubleLazyEdge * > dlazy_edge_vector_;
	mutable InPlaceIntListOP aa_submatrix_history_list_;

	//no default constructor, uncopyable
	DoubleLazyInteractionGraph();
	DoubleLazyInteractionGraph( DoubleLazyInteractionGraph const & );
	DoubleLazyInteractionGraph & operator = ( DoubleLazyInteractionGraph const & );

};

inline
DoubleLazyEdge const * DoubleLazyNode::get_incident_dlazy_edge( int index ) const
{
	return static_cast< DoubleLazyEdge const * >  (get_incident_edge( index ));
}

inline
DoubleLazyEdge * DoubleLazyNode::get_incident_dlazy_edge( int index )
{
	return static_cast< DoubleLazyEdge * >  (get_incident_edge( index ));
}

inline
DoubleLazyNode const * DoubleLazyNode::get_adjacent_dlazy_node( int index ) const
{
	return static_cast< DoubleLazyNode const * > (get_adjacent_node( index ));
}

inline
DoubleLazyNode * DoubleLazyNode::get_adjacent_dlazy_node( int index )
{
	return static_cast< DoubleLazyNode * > (get_adjacent_node( index ));
}

inline
DoubleLazyInteractionGraph const * DoubleLazyNode::get_dlazy_ig_owner() const
{
	return static_cast< DoubleLazyInteractionGraph const * > (get_owner());
}

inline
DoubleLazyInteractionGraph * DoubleLazyNode::get_dlazy_ig_owner()
{
	return static_cast< DoubleLazyInteractionGraph * > (get_owner());
}


inline
DoubleLazyNode const * DoubleLazyEdge::get_dlazy_node( int index ) const
{
	return static_cast< DoubleLazyNode const * > (get_node( index ));
}

inline
DoubleLazyNode * DoubleLazyEdge::get_dlazy_node( int index )
{
	return static_cast< DoubleLazyNode * > (get_node( index ));
}

inline
DoubleLazyInteractionGraph const * DoubleLazyEdge::get_dlazy_ig_owner() const
{
	return static_cast< DoubleLazyInteractionGraph const * > (get_owner());
}

inline
DoubleLazyInteractionGraph * DoubleLazyEdge::get_dlazy_ig_owner()
{
	return static_cast< DoubleLazyInteractionGraph * > (get_owner());
}

inline
void
DoubleLazyEdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state,
	SparseMatrixIndex const & nodes_new_state_sparse_info
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = curr_state_energy;

	get_dlazy_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state,
		nodes_new_state_sparse_info
	);

	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin DoubleLazyNode::project_deltaE_for_substitution
///
/// @brief
///
/// @detailed
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
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
inline
core::PackerEnergy
DoubleLazyNode::project_deltaE_for_substitution
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

	//int aa_neighb_linear_index_offset = aa_offsets_for_edges_.
	//	index(1, 1, alt_state_sparse_mat_info_.get_aa_type() ) - 1;

	//int alt_state_num_states_per_aa_type =
	//	get_num_states_for_aa_type( alt_state_sparse_mat_info_.get_aa_type() );
	//int alt_state_for_aa_type_minus_1 =
	//	alt_state_sparse_mat_info_.get_state_ind_for_this_aa_type() - 1;
	//int nstates_offset =
	//	num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1) - 1;

	for ( int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes(); ++ii ) {

		alternate_state_two_body_energies_[ ii ] = get_incident_dlazy_edge(ii)->
			get_two_body_energy_smi(
				neighbors_curr_state_[ ii ], alternate_state_,
				neighbors_curr_state_sparse_info_[ ii ], alt_state_sparse_mat_info_ );
	}

	for (int ii = get_num_edges_to_smaller_indexed_nodes() + 1;	ii <= get_num_incident_edges(); ++ii) {
		alternate_state_two_body_energies_[ ii ] = get_incident_dlazy_edge(ii)->
			get_two_body_energy_smi(
				alternate_state_, neighbors_curr_state_[ ii ],
				alt_state_sparse_mat_info_, neighbors_curr_state_sparse_info_[ ii ] );
	}

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
	}
	return alternate_state_total_energy_ - curr_state_total_energy_;


}

////////////////////////////////////////////////////////////////////////////////
/// @begin DoubleLazyNode::get_sparse_mat_info_for_curr_state
///
/// @brief
///
/// @detailed
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
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
inline
SparseMatrixIndex const &
DoubleLazyNode::get_sparse_mat_info_for_curr_state() const
{
	return get_sparse_mat_info_for_state( current_state_ );
}


////////////////////////////////////////////////////////////////////////////////
/// @begin DoubleLazyNode::acknowledge_neighbors_state_substitution
///
/// @brief
///
/// @detailed
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
/// @authors apl
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
inline
void
DoubleLazyNode::acknowledge_neighbors_state_substitution(
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


}
}
}

#endif
