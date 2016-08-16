// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/LinearMemoryInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_LinearMemoryInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_LinearMemoryInteractionGraph_hh

// Unit Headers
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>
#include <utility/recent_history_queue.hh>


namespace core {
namespace pack {
namespace interaction_graph {

/// @brief for storing three peices of associated data describing
/// the recent history structure on a LinearMemNode.
struct history_queue_struct
{
	int more_recent_ptr;
	int state_in_rh;
	int more_ancient_ptr;
};

class LinearMemNode : public OnTheFlyNode
{
public:
	/// @brief main constructor, no default ctor, uncopyable
	LinearMemNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states
	);

	/// @brief virtual dstor
	virtual ~LinearMemNode();

	//virtual methods inherited from NodeBase
	/// @brief linmem ig does not have to do anything before sim annealing begins
	virtual void prepare_for_simulated_annealing();
	/// @brief write internal energy and bookkeeping data to standard out
	virtual void print() const;
	/// @brief state 0 represents the unasigned state
	virtual bool state_unassigned() const { return current_state_ == 0;}
	/// @brief return the total energy for this node; includes full energies to neighboring residues
	virtual core::PackerEnergy get_totalE() const { return curr_state_total_energy_;}

	/// @brief set to state 0
	void assign_zero_state();
	/// @brief set to a particular state -- updates the energies internally
	void assign_state(int new_state);
	/// @brief first half of an entire-graph state assignment that avoids unnecessary energy updates
	/// as produced in assign_state.  Adjust all node's states first, then update all energies.
	void partial_assign_state( int new_state );
	/// @brief second half of the entire-graph state assignment.
	void complete_state_assignment();

	/// @brief return the index of the currently assigned state
	inline
	int get_current_state() const
	{ return current_state_; }


	/// @brief return the "recent state id" for the currently assigned state
	inline
	int get_curr_state_recent_state_id() const
	{ return rhq_.head_of_queue(); }

	/// @brief return the one-body energy for the currently assigned state
	inline
	core::PackerEnergy get_one_body_energy_current_state() const
	{ return curr_state_one_body_energy_; }

	/// @brief compute the change in energy induced by substituting the currently assigned state with
	/// some alternative state
	core::PackerEnergy project_deltaE_for_substitution
	(
		int alternate_state,
		core::PackerEnergy & prev_node_energy
	);

	/// @brief proceed to change the currently assigned state to the alternative state considered in
	/// the last call to project_deltaE_for_substitution
	void commit_considered_substitution();

	/// @brief update bookkeeping info to acknolwedge that the last alternative state considered was not
	/// in fact chosen for the currently assigned state.
	void acknowledge_last_substititon_not_committed();

	/// @brief comupute the interaction energy between the currently assigned state on this residue and
	/// the currently assigned state on a neighboring residue, identified by the index of the edge connecting
	/// the two
	core::PackerEnergy compute_pair_energy_for_current_state( int edge_making_energy_request );
	/// @brief comupute the interaction energy between the alternate state being considered on this residue and
	/// the currently assigned state on a neighboring residue, identified by the index of the edge connecting
	/// the two
	core::PackerEnergy compute_pair_energy_for_alternate_state( int edge_making_energy_request );

	inline
	void acknowledge_neighbors_state_substitution(
		int edge_to_altered_neighbor,
		core::PackerEnergy new_edge_energy,
		int other_node_new_state,
		SparseMatrixIndex const & other_node_new_state_sparse_info,
		int other_node_recent_history_index
	);

	void
	acknowledge_neighbors_partial_state_substitution(
		int edge_to_altered_neighbor,
		int other_node_new_state,
		SparseMatrixIndex const & other_node_new_state_sparse_info,
		int other_state_recent_history_index
	);

	inline
	SparseMatrixIndex const &
	get_sparse_mat_info_for_curr_state() const;

	void set_recent_history_size( int num_states_to_maintain_in_recent_history );
	int get_recent_history_size() const;

	void print_internal_energies() const;

	void update_internal_energy_sums();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

protected:

	inline
	conformation::Residue const &
	get_current_rotamer()
	{
		return get_rotamer( current_state_ );
	}

	//Hooks for SASANode< V, E, G > class
	//i.e. for non-PD classes which are templated to use either PD or OTF graphs
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

	void update_internal_vectors();

	inline
	LinearMemEdge const * get_incident_linmem_edge( int index ) const;

	inline
	LinearMemEdge * get_incident_linmem_edge( int index );

	inline
	LinearMemNode const * get_adjacent_linmem_node( int index ) const;

	inline
	LinearMemNode * get_adjacent_linmem_node( int index );

	inline
	LinearMemoryInteractionGraph const *
	get_linmem_ig_owner() const;

	inline
	LinearMemoryInteractionGraph *
	get_linmem_ig_owner();

	int update_recent_history( int state );

private:
	/// Data
	utility::recent_history_queue rhq_;

	ObjexxFCL::FArray3D< unsigned char > aa_neighbors_for_edges_;
	utility::vector1< int > neighbors_curr_state_;
	utility::vector1< int > neighbors_state_recent_history_index_;
	utility::vector1< SparseMatrixIndex > neighbors_curr_state_sparse_info_;

	int current_state_;
	SparseMatrixIndex curr_state_sparse_mat_info_;
	core::PackerEnergy curr_state_one_body_energy_;
	core::PackerEnergy curr_state_total_energy_;
	utility::vector1< core::PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
	SparseMatrixIndex alt_state_sparse_mat_info_;
	core::PackerEnergy alternate_state_one_body_energy_;
	core::PackerEnergy alternate_state_total_energy_;
	utility::vector1< core::PackerEnergy > alternate_state_two_body_energies_;

	bool alternate_state_is_being_considered_;
	bool already_prepped_for_simA_;

	ObjexxFCL::FArray1D_int accepted_rejected_substitution_history_;
	int accepted_history_head_;
	int num_recently_accepted_;
	bool filled_substitution_history_;

	static const int ACCEPTED = 1;
	static const int REJECTED = 0;

	static const int ACCEPTANCE_REJECTION_HISTORY_LENGTH = 100;
	static const int THRESHOLD_ACCEPTANCE_RATE_FOR_RPE_STORAGE = 10;
};

class LinearMemEdge : public OnTheFlyEdge
{
public:
	LinearMemEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind
	);

	virtual ~LinearMemEdge();

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
	virtual void prepare_for_simulated_annealing();

	core::PackerEnergy get_current_two_body_energy() const;

	void acknowledge_state_change(
		int node_ind,
		int new_state,
		SparseMatrixIndex const & new_state_sparse_info,
		int bumped_recent_history_index,
		int new_state_recent_history_index,
		core::PackerEnergy & new_energy
	);
	void acknowledge_state_zeroed( int node_ind );

	void acknowledge_partial_state_change(
		int node_ind,
		int new_state,
		SparseMatrixIndex const & new_state_sparse_info,
		int bumped_recent_history_index,
		int new_state_recent_history_index
	);

	core::PackerEnergy get_energy_following_partial_state_assignment();

	void reset_state_energies(
		int node_index,
		int state,
		int recent_history_id
	);

	core::PackerEnergy
	get_energy_for_alt_state
	(
		bool store_rpes,
		int changing_node_index,
		int alternate_state,
		int alternate_state_recent_history_index,
		int other_node_curr_state,
		int other_node_state_recent_history_index
	);

	inline
	void acknowledge_substitution(
		int substituted_node_index,
		core::PackerEnergy const curr_state_energy,
		int nodes_new_state,
		SparseMatrixIndex const & nodes_new_state_sparse_info,
		int bumped_recent_history_index,
		int new_state_recent_history_index,
		int neighbors_curr_state
	);

	int get_two_body_table_size() const;
	virtual void declare_energies_final();

	void print_current_energy() const;


	ObjexxFCL::FArray2D< unsigned char > const & get_sparse_aa_neighbor_info( );
	//ObjexxFCL::FArray1D_int const & get_second_node_num_states_per_aa();

	static core::PackerEnergy const NOT_YET_COMPUTED_ENERGY; //an energy lower than any RPE could ever be

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	virtual void set_edge_weight( Real weight );

protected:

	//Hooks for SASAEdge< V, E, G > class
	void declare_energies_final_no_deletion() { LinearMemEdge::prepare_for_simulated_annealing(); }
	void prepare_for_simulated_annealing_no_deletion() { LinearMemEdge::prepare_for_simulated_annealing(); }
	bool pd_edge_table_all_zeros() const { return false;}


private:

	inline LinearMemNode const * get_linmem_node( int index ) const;
	inline LinearMemNode * get_linmem_node( int index );
	inline LinearMemoryInteractionGraph const * get_linmem_ig_owner() const;
	inline LinearMemoryInteractionGraph * get_linmem_ig_owner();

	void
	handle_bumped_recent_history_state_for_node(
		int node_substituted,
		int node_not_substituted,
		int bumped_recent_history_index
	);

	void store_curr_state_energy();

	void wipe( int node );

	bool store_rpes_[ 2 ];
	ObjexxFCL::FArray2D< core::PackerEnergy > stored_rpes_[ 2 ];
	ObjexxFCL::FArray2D< unsigned char > sparse_aa_neighbors_;
	core::PackerEnergy curr_state_energy_;
	core::PackerEnergy alt_state_energy_;
	bool partial_state_assignment_;
	bool preped_for_sim_annealing_;

	//no default constructor, uncopyable
	LinearMemEdge();
	LinearMemEdge( LinearMemEdge const & );
	LinearMemEdge & operator = ( LinearMemEdge const & );

};

class LinearMemoryInteractionGraph : public OnTheFlyInteractionGraph
{
public:
	LinearMemoryInteractionGraph( int numNodes );
	virtual ~LinearMemoryInteractionGraph();

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
	virtual core::PackerEnergy commit_considered_substitution();
	virtual core::PackerEnergy get_energy_current_state_assignment();
	virtual int get_edge_memory_usage() const;
	virtual void print_current_state_assignment() const;
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE );
	virtual core::PackerEnergy get_energy_sum_for_vertex_group( int group_id );
	virtual void prepare_for_simulated_annealing();

	void set_recent_history_size( Size recent_history_size );
	Size get_recent_history_size() const;

	//bool build_sc_only_rotamer() const;

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

protected:

	virtual NodeBase* create_new_node( int node_index, int num_states);
	virtual EdgeBase* create_new_edge( int index1, int index2);

	//Hooks for SASAInterationGraph< V, E, G >
	core::PackerEnergy get_energy_PD_current_state_assignment();
	void update_internal_energy_totals();

	inline
	LinearMemNode const * get_linmem_node(int index) const
	{
		return static_cast< LinearMemNode const * > (get_node( index ));
	}

	inline
	LinearMemNode * get_linmem_node(int index)
	{
		return static_cast< LinearMemNode * > (get_node( index ));
	}

private:
	/// @brief Set the recent history size for all nodes in the graph
	void set_recent_history_sizes();

	bool first_time_prepping_for_simA_;
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;
	int recent_history_size_;

	bool have_not_committed_last_substitution_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	LinearMemoryInteractionGraph();
	LinearMemoryInteractionGraph( LinearMemoryInteractionGraph const & );
	LinearMemoryInteractionGraph & operator = ( LinearMemoryInteractionGraph const & );
};


inline
LinearMemoryInteractionGraph const *
LinearMemNode::get_linmem_ig_owner() const
{
	return static_cast< LinearMemoryInteractionGraph const * > (get_owner());
}

inline
LinearMemoryInteractionGraph *
LinearMemNode::get_linmem_ig_owner()
{
	return static_cast< LinearMemoryInteractionGraph * > (get_owner());
}

inline
LinearMemEdge const *
LinearMemNode::get_incident_linmem_edge( int index ) const
{
	return static_cast< LinearMemEdge const * > (get_incident_edge( index ));
}

inline
LinearMemEdge *
LinearMemNode::get_incident_linmem_edge( int index )
{
	return static_cast< LinearMemEdge * > (get_incident_edge( index ));
}

inline
LinearMemNode const *
LinearMemNode::get_adjacent_linmem_node( int index ) const
{
	return static_cast< LinearMemNode const * > (get_adjacent_node( index ));
}

inline
LinearMemNode *
LinearMemNode::get_adjacent_linmem_node( int index )
{
	return static_cast< LinearMemNode * > (get_adjacent_node( index ));
}

inline
LinearMemNode const *
LinearMemEdge::get_linmem_node( int index ) const
{
	return static_cast< LinearMemNode const * > (get_node( index ));
}

inline
LinearMemNode *
LinearMemEdge::get_linmem_node( int index )
{
	return static_cast< LinearMemNode * > (get_node( index ));
}


inline
LinearMemoryInteractionGraph const *
LinearMemEdge::get_linmem_ig_owner() const
{
	return static_cast< LinearMemoryInteractionGraph const * > (get_owner());
}

inline
LinearMemoryInteractionGraph *
LinearMemEdge::get_linmem_ig_owner()
{
	return static_cast< LinearMemoryInteractionGraph * > (get_owner());
}


inline
SparseMatrixIndex const &
LinearMemNode::get_sparse_mat_info_for_curr_state() const
{
	return get_sparse_mat_info_for_state( current_state_ );
}

inline
void
LinearMemEdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state,
	SparseMatrixIndex const & nodes_new_state_sparse_info,
	int bumped_recent_history_index,
	int new_state_recent_history_index,
	int neighbors_curr_state
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	handle_bumped_recent_history_state_for_node(
		node_substituted,
		node_not_substituted,
		bumped_recent_history_index );

	curr_state_energy_ = curr_state_energy;
	if ( neighbors_curr_state != 0 ) {
		stored_rpes_[ node_substituted ]
			( neighbors_curr_state, new_state_recent_history_index ) =
			curr_state_energy_;
	}

	get_linmem_node( node_not_substituted )->
		acknowledge_neighbors_state_substitution
		(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state,
		nodes_new_state_sparse_info,
		new_state_recent_history_index
	);

	return;
}

inline
void
LinearMemNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	core::PackerEnergy new_edge_energy,
	int other_node_new_state,
	SparseMatrixIndex const & other_node_new_state_sparse_info,
	int other_node_recent_history_index
)
{
	curr_state_total_energy_ +=
		new_edge_energy - curr_state_two_body_energies_[ edge_to_altered_neighbor ];
	curr_state_two_body_energies_[ edge_to_altered_neighbor ] = new_edge_energy;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	neighbors_curr_state_sparse_info_[ edge_to_altered_neighbor ]  =
		other_node_new_state_sparse_info;
	neighbors_state_recent_history_index_[ edge_to_altered_neighbor ] =
		other_node_recent_history_index;
	return;
}


}
}
}

#endif

