// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_interaction_graph_FASTERInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_FASTERInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/FASTERInteractionGraph.fwd.hh>

// Package Headers

#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>

//STL Headers
#include <list>

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class FASTERNode;
class FASTEREdge;

class FASTERNode : public PrecomputedPairEnergiesNode
{
public:
	typedef PrecomputedPairEnergiesNode parent;
public:
	FASTERNode(InteractionGraphBase * owner, int node_id, int num_states);
	virtual ~FASTERNode();
	virtual void print() const;

	void update_one_body_energy( int state, core::PackerEnergy energy );
	virtual void update_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies );
	void add_to_one_body_energy( int state, core::PackerEnergy energy );
	virtual void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies );
	virtual void zero_one_body_energies();
	core::PackerEnergy get_one_body_energy( int state );

	virtual void prepare_for_simulated_annealing();
	void prepare_for_FASTER();
	//virtual unsigned int getMemoryUsageInBytes() const;

	void assign_zero_state();
	virtual bool state_unassigned() const { return current_state_ == 0; }
	void assign_state(int new_state);
	int get_current_state() const;
	core::PackerEnergy get_one_body_energy_current_state() const;
	core::PackerEnergy project_deltaE_for_substitution
	(
		int alternate_state,
		core::PackerEnergy & prev_node_energy
	);
	void commit_considered_substitution();

	inline
	void acknowledge_neighbors_state_substitution(
		int edge_to_altered_neighbor,
		core::PackerEnergy new_edge_energy,
		int other_node_new_state);

	void print_internal_energies() const;
	void update_internal_energy_sums();

	/*
	void prepare_to_write_to_file();
	void initialize_aa_for_state_array();
	void clean_up_after_writing_to_file();
	void prepare_to_read_energies_from_file( int num_states_for_node_in_file );
	void clean_up_after_reading_energies_from_file();

	void set_aa_for_file_state(int file_state, int aa );
	void set_instance_state_correspondence( int instance_state, int state_from_file );
	int get_correspondence_for_state( int instance_state );

	int get_num_rots_absent_from_file();
	void get_absent_rots( ObjexxFCL::FArray1DB_int & rots_absent );

	int get_num_states_in_file();
	int & get_aatypes_for_file_states();

	int & get_aatypes_for_states();
	int & get_num_file_states_for_aa();
	int & get_file_states_2_instance_states_array();

	bool get_node_corresponded_to_file_node();
	 */

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

protected:
	void update_internal_vectors();

	inline FASTEREdge const * get_incident_faster_edge( int index ) const;
	inline FASTEREdge       * get_incident_faster_edge( int index );

	inline FASTERNode const * get_adjacent_faster_node( int index ) const;
	inline FASTERNode       * get_adjacent_faster_node( int index );

	inline FASTERInteractionGraph const * get_fasterig_owner() const;
	inline FASTERInteractionGraph       * get_fasterig_owner();

public:

	// Methods for FASTER
	void get_total_energy_in_curr_state_assignment_for_all_states();
	void partial_assign_state_with_lowest_one_body_energy();
	void partial_assign_relaxed_state( Real probability );
	void partial_assign_state( int state );
	void acknowledge_neighbors_partial_state_assignment(
		int which_neighbor,
		int neighbors_new_state
	);
	void complete_partial_state_assignment();
	void acknowledge_neighbors_perturbed_state(
		int which_neighbor,
		int const neighbors_perturbed_state
	);

	void prepare_for_perturbation();

	void set_no_longer_perturbed();

	void set_perturbed_state( int perturbed_state );


	void relax_neighbors();

	void relax();

	void reset_relaxed_for_neighbors();
	void reset_relaxed();
	void tell_neighbors_to_prep_for_relaxation();
	void prep_for_neighbors_perturbation();
	void relax_after_neighbors_perturbation();
	core::PackerEnergy get_deltaE_for_relaxed_state_following_perturbation();
	core::PackerEnergy get_one_body_energy_for_relaxed_state() const;
	int get_relaxed_state() const;
	int get_random_neighbor();

private:
	// Data

	std::vector< core::PackerEnergy > one_body_energies_;

	std::vector< int > neighbors_curr_state_;
	std::vector< ObjexxFCL::FArray2A< core::PackerEnergy > > edge_matrix_ptrs_;

	int current_state_;
	core::PackerEnergy curr_state_one_body_energy_;
	core::PackerEnergy curr_state_total_energy_;
	std::vector< core::PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
	core::PackerEnergy alternate_state_one_body_energy_;
	core::PackerEnergy alternate_state_total_energy_;
	std::vector< core::PackerEnergy > alternate_state_two_body_energies_;

	bool alternate_state_is_being_considered_;

	//variables for FASTER
	bool have_prepared_once_for_FASTER_;
	bool in_FASTER_mode_;
	bool perturbed_;
	bool have_relaxed_since_neighbors_perturbation_;
	bool have_contributed_deltaE_following_perturbation_;
	int relaxed_state_;
	std::vector< core::PackerEnergy > state_energies_in_current_state_assignment_;
	std::vector< core::PackerEnergy > state_energies_in_current_context_;
	std::vector< bool > neighbor_relaxed_in_sBR_;
	std::vector< core::PackerEnergy > perturbed_two_body_energies_;

	//no default constructor, uncopyable
	FASTERNode();
	FASTERNode( FASTERNode const & );
	FASTERNode & operator = ( FASTERNode const & );
};

class FASTEREdge : public PrecomputedPairEnergiesEdge
{
public:
	typedef PrecomputedPairEnergiesEdge parent;

public:
	FASTEREdge(InteractionGraphBase* owner, int first_node_ind, int second_node_ind);
	virtual ~FASTEREdge();
	virtual void set_sparse_aa_info(ObjexxFCL::FArray2_bool const & ) {}
	virtual bool get_sparse_aa_info( int, int) const {return true;} //"all amino acids are neighbors"
	virtual void add_to_two_body_energy(int const, int const, core::PackerEnergy const);
	virtual void add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array );
	virtual
	void set_two_body_energy(int const, int const, core::PackerEnergy const);
	virtual
	void clear_two_body_energy(int const, int const);
	virtual core::PackerEnergy get_two_body_energy( int const, int const ) const;

	virtual void force_aa_neighbors(int, int) {} //all aa's are already neighbors -- dense representation
	virtual void force_all_aa_neighbors() {} //same thing


	virtual void declare_energies_final();
	virtual void prepare_for_simulated_annealing();
	virtual void prepare_for_FASTER();
	//virtual unsigned int getMemoryUsageInBytes() const;

	core::PackerEnergy get_current_two_body_energy();

	void acknowledge_partial_state_assignment(
		int node,
		int new_state
	);

	void acknowledge_state_change(
			int node_ind,
			int new_state,
			core::PackerEnergy & new_energy
	);
	void acknowledge_state_zeroed( int node_ind );

	static
	inline
	core::PackerEnergy get_alternate_state_energy(
		int first_node_state,
		int second_node_state,
		ObjexxFCL::FArray2< core::PackerEnergy > & edge_energy_table
	);

	inline void acknowledge_substitution(
		int substituted_node_index,
		core::PackerEnergy const curr_state_energy,
		int nodes_new_state
	);

	ObjexxFCL::FArray2A< core::PackerEnergy > get_edge_table_ptr();

	int get_two_body_table_size() const;

	//ObjexxFCL::FArray2Da< core::PackerEnergy > get_lonode_edge_table_ptr();
	//ObjexxFCL::FArray2Da< core::PackerEnergy > get_hinode_edge_table_ptr();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	//methods for FASTER
	core::PackerEnergy get_two_body_energies_for_relaxed_states();

	core::PackerEnergy get_curr_state_energy_following_partial_state_assignment();

	core::PackerEnergy
	get_energy_for_perturbed_state(
		int node,
		int nodes_perturbed_state
	);

	void
	acknowledge_perturbed_state(
		int node,
		int neighbors_context_state
	);

	void acknowledge_participation_in_perturbation();
	core::PackerEnergy get_deltaE_for_perturbation();
	core::PackerEnergy get_deltaE_for_neighbor_following_perturbation( int node_index );

	void set_edge_weight( Real weight );

	void
	swap_edge_energies(
		ObjexxFCL::FArray2D< core::PackerEnergy > & new_edge_table
	);

private:
	inline
	FASTERNode* get_faster_node( int index ) const
	{
		return (FASTERNode*) get_node( index ); //c-style cast since static_cast won't compile
	}

	inline
	FASTERInteractionGraph* get_fasterig_owner() const
	{
		return (FASTERInteractionGraph*) get_owner(); //c-style cast since static_cast won't compile
	}

	//ObjexxFCL::FArray2D< core::PackerEnergy > lonode_two_body_energies_; //Dense matrix
	//ObjexxFCL::FArray2D< core::PackerEnergy > hinode_two_body_energies_; //Dense matrix

private:
	// DATA
	ObjexxFCL::FArray2D< core::PackerEnergy > two_body_energies_;
	core::PackerEnergy curr_state_energy_;
	bool energies_updated_since_last_prep_for_simA_;
	bool have_contributed_deltaE_following_perturbation_;
	bool partial_state_assignment_;

private:
	// Unimplemented functions
	//no default constructor, uncopyable
	FASTEREdge();
	FASTEREdge( FASTEREdge const & );
	FASTEREdge & operator = ( FASTEREdge const & );
};

class FASTERInteractionGraph : public PrecomputedPairEnergiesInteractionGraph
{
public:
	typedef PrecomputedPairEnergiesInteractionGraph parent;

public:
	FASTERInteractionGraph(int num_nodes);
	virtual void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

	virtual core::PackerEnergy get_one_body_energy_for_node_state( int node, int state);

	//virtual void set_num_aatypes(int) {}
	virtual int  get_num_aatypes() const {return 1;}

	virtual void blanket_assign_state_0();
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state);
	int get_current_state_for_node( int node_ind ) const;
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states);
	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	);

	/// @brief Accepts (commits) the state change previously considered in a call to
	/// consider_substitution and returns the energy of the entire graph
	virtual core::PackerEnergy commit_considered_substitution();

	/// @brief removes all accumulated numerical drift and returns the
	/// energy for the current state assignment.
	virtual core::PackerEnergy get_energy_current_state_assignment();
	/// @brief returns the number of core::PackerEnergys used in all edge two-body energy tables
	virtual int get_edge_memory_usage() const;

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	//methods for FASTER optimization
	virtual void prepare_for_FASTER();
	void assign_BMEC();
	void relax_in_current_context();
	core::PackerEnergy get_energy_following_relaxation();
	void commit_relaxation();
	void probabilistically_commit_relaxation( Real probability );
	void get_current_network_state( ObjexxFCL::FArray1_int & netstate );
	virtual core::PackerEnergy perturb_sBR_and_relax( int node, int perturbed_state );
	virtual core::PackerEnergy perturb_dBR_and_relax( int node1, int perturbed_state1, int node2, int perturbed_state2 );
	int get_random_neighbor_for_node( int node );
	void reject_perturbation();

	/// @brief outputs the current state for each node, useful for debugging
	virtual void print_current_state_assignment() const;
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy ) {};

	/// @brief a user may define subsets of the vertex set for which they would like to
	/// know the internal energy sum.
	virtual core::PackerEnergy get_energy_sum_for_vertex_group( int group_id );

	/// @brief Swap the contents of the edge pair-energy table with the contents in the new_edge_table
	void
	swap_edge_energies( int node1, int node2, ObjexxFCL::FArray2D< core::PackerEnergy > & new_edge_table );

protected:
	//virtual unsigned int getMemoryUsageInBytes() const;

	virtual NodeBase* create_new_node( int node_index, int num_states);
	virtual EdgeBase* create_new_edge( int index1, int index2);

	/// @brief removes numerical drift that can accumulate over the course of
	/// many state assignment changes within simulated annealing
	void update_internal_energy_totals();

	inline
	FASTERNode const * get_faster_node(int index) const
	{
		return static_cast< FASTERNode const * > (get_node( index ));
	}

	inline
	FASTERNode       * get_faster_node(int index)
	{
		return static_cast< FASTERNode       * > (get_node( index ));
	}

	inline
	FASTEREdge const * get_faster_edge(int node1, int node2) const
	{
		return static_cast< FASTEREdge const * > (find_edge(node1, node2));
	}

	inline
	FASTEREdge       * get_faster_edge(int node1, int node2)
	{
		return static_cast< FASTEREdge       * > (find_edge(node1, node2));
	}


private:
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;

	//variables for FASTER
	bool sBR_;
	bool dBR_;
	int relaxed1_;
	int relaxed2_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	FASTERInteractionGraph();
	FASTERInteractionGraph( FASTERInteractionGraph const & );
	FASTERInteractionGraph & operator = ( FASTERInteractionGraph const & );
};

inline
FASTEREdge const * FASTERNode::get_incident_faster_edge( int index ) const
{
	return static_cast< FASTEREdge const * > (get_incident_edge( index ));
}

inline
FASTEREdge       * FASTERNode::get_incident_faster_edge( int index )
{
	return static_cast< FASTEREdge       * > (get_incident_edge( index ));
}

inline
FASTERNode const * FASTERNode::get_adjacent_faster_node( int index ) const
{
	return static_cast< FASTERNode const * > (get_adjacent_node( index ));
}

inline
FASTERNode       * FASTERNode::get_adjacent_faster_node( int index )
{
	return static_cast< FASTERNode * > (get_adjacent_node( index ));
}

inline
FASTERInteractionGraph const * FASTERNode::get_fasterig_owner() const
{
	return static_cast< FASTERInteractionGraph const * > (get_owner());
}

inline
FASTERInteractionGraph       * FASTERNode::get_fasterig_owner()
{
	return static_cast< FASTERInteractionGraph       * > (get_owner());
}

/// @brief updates bookkeeping arrays for when a neighbor has changed its state
///
/// @param edge_to_altered_neighbor - [in] - the index for the edge that connects
/// 	this node to the node that just changed its state
/// @param new_edge_energ - [in] - the pair energy between this node in its current
///	state and the new state of the node that just changed its state
/// @param other_node_new_state - [in] - the state the neighbor just adopted
inline
void FASTERNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	core::PackerEnergy new_edge_energy,
	int other_node_new_state
)
{

	curr_state_total_energy_ += new_edge_energy - curr_state_two_body_energies_[edge_to_altered_neighbor];
	curr_state_two_body_energies_[edge_to_altered_neighbor] = new_edge_energy;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	//neighbors_curr_state_plus_offset_[ edge_to_altered_neighbor ]
	//	= other_node_new_state + neighbors_rotindex_offset_[ edge_to_altered_neighbor ];
	return;
}

/// @brief static method that looks up the two body energy when the
/// node with the smaller index on an edge is considering an alternate state
///
/// @param first_node_alt_state - [in] - the alternate state for the lower-indexed node
/// @param second_node_orig_state - [in] - the current state for the higher-indexed node
/// @param edge_energy_table - [in] - the proxy FArray pointing at the edge table
/// 	connecting the two nodes.
inline
core::PackerEnergy
FASTEREdge::get_alternate_state_energy(
	int first_node_state,
	int second_node_state,
	ObjexxFCL::FArray2< core::PackerEnergy > & edge_energy_table
)
{
	if (first_node_state == 0 || second_node_state == 0) {
		return core::PackerEnergy( 0.0 );
	} else {
		return edge_energy_table( second_node_state, first_node_state );
	}
}

/// @brief update bookkeeping information when one of the nodes an edge is incident
/// upon changes state
///
/// @param substituted_node_index - [in] - index of the node that chagned its state
/// @param curr_state_energy - [in] - the two body energy given the new state
/// @param nodes_new_state - [in] - the state the node just transitioned into
/// @param nodes_new_state_sparse_info - [in] - sparse matrix info for the new state
inline
void
FASTEREdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = curr_state_energy;

	get_faster_node( node_not_substituted )->acknowledge_neighbors_state_substitution(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state
	);
}

/// @brief returns the change in energy that would be induced by switching this node
/// from its current state into another state
///
/// iterates across the incident edges for a node in two phases:
/// in the first phase, it examines edges leading to higher-indexed nodes
/// in the second phase, it examines edges leading to smaller-indexed nodes.
/// for cache efficiency, all of the amino-acid-neighbor-offset information
/// that each edge calculates is stored on the nodes themselves.  The edges
/// are never touched; rather, their private information is stored on the nodes
/// and handed to static member functions of the DensePDEdge class.  This "store
/// edge information on the nodes" strategy gives me performance equivalent
/// to the previous energy2b lookup tables.
///
/// @param alternate_state - [in] - the alternate state to consider
/// @param previous_energy_for_node - [out] - the old energy1b/energy2b sum for this
/// node; used by simulate annealing.
///
inline
core::PackerEnergy
FASTERNode::project_deltaE_for_substitution(
	int alternate_state,
	core::PackerEnergy & prev_energy_for_node
)
{

	alternate_state_is_being_considered_ = true;
	//std::cerr << "proj_deltaE: node -  " << get_node_index()
	// << " alt state " << alternate_state << "...";

	alternate_state_ = alternate_state;
	alternate_state_one_body_energy_ = one_body_energies_[ alternate_state ];
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
	prev_energy_for_node = curr_state_total_energy_;

	for ( int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes(); ++ii ) {
		alternate_state_two_body_energies_[ ii ] =
			FASTEREdge::get_alternate_state_energy(
				neighbors_curr_state_[ii],
				alternate_state_,
				edge_matrix_ptrs_[ii]
		);
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
	}

	for ( int ii = get_num_edges_to_smaller_indexed_nodes() + 1; ii <= get_num_incident_edges(); ++ii ) {
		alternate_state_two_body_energies_[ ii ] =
			FASTEREdge::get_alternate_state_energy(
				alternate_state_,
				neighbors_curr_state_[ii],
				edge_matrix_ptrs_[ii]
		);
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
	}

	//std::cerr<< "..done" << std::endl;

	return alternate_state_total_energy_ - curr_state_total_energy_;

}


} // namespace interaction_graph
} // namespace pack
} // namespace core

#endif
