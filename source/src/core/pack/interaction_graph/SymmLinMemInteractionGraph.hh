// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/SymmLinMemInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_SymmLinMemInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_SymmLinMemInteractionGraph_hh

// Unit Headers
#include <core/pack/interaction_graph/SymmLinMemInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>
#include <utility/recent_history_queue.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class SymmLinearMemNode : public SymmOnTheFlyNode
{
public:
	/// @brief main constructor, no default ctor, uncopyable
	SymmLinearMemNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states
	);

	/// @brief virtual dstor
	virtual ~SymmLinearMemNode();

	//virtual methods inherited from NodeBase
	/// @brief symmlinmem ig does not have to do anything before sim annealing begins
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
		int other_node_recent_history_index
	);

	void
	acknowledge_neighbors_partial_state_substitution(
		int edge_to_altered_neighbor,
		int other_node_new_state,
		int other_state_recent_history_index
	);

	void set_recent_history_size( int num_states_to_maintain_in_recent_history );
	int get_recent_history_size() const;

	void print_internal_energies() const;

	void update_internal_energy_sums();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

private:

	void update_internal_vectors();

public:
	inline
	SymmLinearMemEdge const * get_incident_symmlinmem_edge( int index ) const;

	inline
	SymmLinearMemEdge * get_incident_symmlinmem_edge( int index );

	inline
	SymmLinearMemNode const * get_adjacent_symmlinmem_node( int index ) const;

	inline
	SymmLinearMemNode * get_adjacent_symmlinmem_node( int index );

	inline
	SymmLinearMemoryInteractionGraph const *
	get_symmlinmem_ig_owner() const;

	inline
	SymmLinearMemoryInteractionGraph *
	get_symmlinmem_ig_owner();

private:

	int update_recent_history( int state );

private:
	/// Data
	utility::recent_history_queue rhq_;

	utility::vector1< int > neighbors_curr_state_;
	utility::vector1< int > neighbors_state_recent_history_index_;

	int current_state_;
	core::PackerEnergy curr_state_one_body_energy_;
	core::PackerEnergy curr_state_total_energy_;
	utility::vector1< core::PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
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

class SymmLinearMemEdge : public SymmOnTheFlyEdge
{
public:
	SymmLinearMemEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind
	);

	virtual ~SymmLinearMemEdge();

	virtual core::PackerEnergy get_two_body_energy( int const node1state, int const node2state) const;

	//virtual methods inherited from EdgeBase
	virtual void prepare_for_simulated_annealing();

	core::PackerEnergy get_current_two_body_energy() const;

	void acknowledge_state_change(
		int node_ind,
		int new_state,
		int bumped_recent_history_index,
		int new_state_recent_history_index,
		core::PackerEnergy & new_energy
	);
	void acknowledge_state_zeroed( int node_ind );

	void acknowledge_partial_state_change(
		int node_ind,
		int new_state,
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
	get_energy_for_alt_state(
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
		int bumped_recent_history_index,
		int new_state_recent_history_index,
		int neighbors_curr_state
	);

	int get_two_body_table_size() const;
	virtual void declare_energies_final();

	void print_current_energy() const;

	static core::PackerEnergy const NOT_YET_COMPUTED_ENERGY; //an energy lower than any RPE could ever be

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	virtual void set_edge_weight( Real weight );

public:

	inline SymmLinearMemNode const * get_symmlinmem_node( int index ) const;
	inline SymmLinearMemNode * get_symmlinmem_node( int index );
	inline SymmLinearMemoryInteractionGraph const * get_symmlinmem_ig_owner() const;
	inline SymmLinearMemoryInteractionGraph * get_symmlinmem_ig_owner();

private:

	void
	handle_bumped_recent_history_state_for_node(
		int node_substituted,
		int node_not_substituted,
		int bumped_recent_history_index
	);

	void store_curr_state_energy();

	void wipe( int node );

private:
	bool store_rpes_[ 2 ];
	ObjexxFCL::FArray2D< core::PackerEnergy > stored_rpes_[ 2 ];
	core::PackerEnergy curr_state_energy_;
	core::PackerEnergy alt_state_energy_;
	bool partial_state_assignment_;
	bool preped_for_sim_annealing_;

	//no default constructor, uncopyable
	SymmLinearMemEdge();
	SymmLinearMemEdge( SymmLinearMemEdge const & );
	SymmLinearMemEdge & operator = ( SymmLinearMemEdge const & );

};

class SymmLinearMemoryInteractionGraph : public SymmOnTheFlyInteractionGraph
{
public:
	SymmLinearMemoryInteractionGraph( int numNodes );
	virtual ~SymmLinearMemoryInteractionGraph();

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
	SymmLinearMemNode const * get_symmlinmem_node(int index) const
	{
		return static_cast< SymmLinearMemNode const * > (get_node( index ));
	}

	inline
	SymmLinearMemNode * get_symmlinmem_node(int index)
	{
		return static_cast< SymmLinearMemNode * > (get_node( index ));
	}

private:
	void set_recent_history_sizes();
	//static int get_cmdline_history_size();

	bool first_time_prepping_for_simA_;
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;
	int recent_history_size_;

	bool have_not_committed_last_substitution_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	SymmLinearMemoryInteractionGraph();
	SymmLinearMemoryInteractionGraph( SymmLinearMemoryInteractionGraph const & );
	SymmLinearMemoryInteractionGraph & operator = ( SymmLinearMemoryInteractionGraph const & );
};


inline
SymmLinearMemoryInteractionGraph const *
SymmLinearMemNode::get_symmlinmem_ig_owner() const
{
	return static_cast< SymmLinearMemoryInteractionGraph const * > (get_owner());
}

inline
SymmLinearMemoryInteractionGraph *
SymmLinearMemNode::get_symmlinmem_ig_owner()
{
	return static_cast< SymmLinearMemoryInteractionGraph * > (get_owner());
}

inline
SymmLinearMemEdge const *
SymmLinearMemNode::get_incident_symmlinmem_edge( int index ) const
{
	return static_cast< SymmLinearMemEdge const * > (get_incident_edge( index ));
}

inline
SymmLinearMemEdge *
SymmLinearMemNode::get_incident_symmlinmem_edge( int index )
{
	return static_cast< SymmLinearMemEdge * > (get_incident_edge( index ));
}

inline
SymmLinearMemNode const *
SymmLinearMemNode::get_adjacent_symmlinmem_node( int index ) const
{
	return static_cast< SymmLinearMemNode const * > (get_adjacent_node( index ));
}

inline
SymmLinearMemNode *
SymmLinearMemNode::get_adjacent_symmlinmem_node( int index )
{
	return static_cast< SymmLinearMemNode * > (get_adjacent_node( index ));
}

inline
SymmLinearMemNode const *
SymmLinearMemEdge::get_symmlinmem_node( int index ) const
{
	return static_cast< SymmLinearMemNode const * > (get_node( index ));
}

inline
SymmLinearMemNode *
SymmLinearMemEdge::get_symmlinmem_node( int index )
{
	return static_cast< SymmLinearMemNode * > (get_node( index ));
}


inline
SymmLinearMemoryInteractionGraph const *
SymmLinearMemEdge::get_symmlinmem_ig_owner() const
{
	return static_cast< SymmLinearMemoryInteractionGraph const * > (get_owner());
}

inline
SymmLinearMemoryInteractionGraph *
SymmLinearMemEdge::get_symmlinmem_ig_owner()
{
	return static_cast< SymmLinearMemoryInteractionGraph * > (get_owner());
}


inline
void
SymmLinearMemEdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state,
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

	get_symmlinmem_node( node_not_substituted )->
		acknowledge_neighbors_state_substitution
		(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state,
		new_state_recent_history_index
	);

	return;
}

inline
void
SymmLinearMemNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	core::PackerEnergy new_edge_energy,
	int other_node_new_state,
	int other_node_recent_history_index
)
{
	curr_state_total_energy_ +=
		new_edge_energy - curr_state_two_body_energies_[ edge_to_altered_neighbor ];
	curr_state_two_body_energies_[ edge_to_altered_neighbor ] = new_edge_energy;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	neighbors_state_recent_history_index_[ edge_to_altered_neighbor ] =
		other_node_recent_history_index;
	return;
}


}
}
}

#endif

