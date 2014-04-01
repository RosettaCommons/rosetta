// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/PDInteractionGraph.hh
/// @brief  Pairwise Decomposable interaction graph class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_PDInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_PDInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/PDInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/SparseMatrixIndex.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>

#include <core/types.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility Headers
//#include <utility/io/all.fwd.hh>

// C++ Headers
// AUTO-REMOVED #include <vector>
#include <list>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1A.hh>



namespace core {
namespace pack {
namespace interaction_graph {

class PDNode;
class PDEdge;
class PDInteractionGraph;

class PDNode : public PrecomputedPairEnergiesNode
{
public:
	PDNode(InteractionGraphBase * owner, int node_id, int num_states);
	virtual ~PDNode();
	/// @brief prints a description of the node and all of it's one-body energies
	virtual void print() const;
	/// @brief sets the amino acid type for each state
	virtual void set_amino_acid_types( std::vector< int > const & );
	/// @brief return the amino acid type for a particular state -- this indexing is of course completely
	/// arbitrary
	virtual
	int aatype_for_state( int state ) const;
	/// @brief returns an FArray & with the number of states for each amino acid type
	utility::vector1< int > const & get_num_states_for_aa_types() const;
	/// @brief update energy to the one-body energy for state
	virtual void update_one_body_energy( int state, core::PackerEnergy energy );
	/// @brief set all the one-body energies for this node
	virtual void update_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies );
	/// @brief adds energy to the one-body energy for state state
	virtual void add_to_one_body_energy( int state, core::PackerEnergy energy );
	/// @brief adds all the energies in energies to the one-body energies for this node
	virtual void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies );
	virtual void zero_one_body_energies();
	/// @brief returns the one body energy for a state
	core::PackerEnergy get_one_body_energy( int state );
	/// @brief prepares node for simulated annealing
	virtual void prepare_for_simulated_annealing();

	/// @brief assigns node's state to it's zero, or "unassigned" state.
	void assign_zero_state();
	virtual bool state_unassigned() const { return current_state_ == 0;}
	/// @brief assigns node a new_state
	void assign_state(int new_state);
	/// @brief returns the state the node is currently assigned
	int get_current_state() const;
	/// @brief returns the one body energy for the state the node is currently assigned
	core::PackerEnergy get_one_body_energy_current_state() const;

	/// @brief
	/// returns the change in energy that would be induced by switching this node
	/// from its current state into another state
	///
	inline core::PackerEnergy project_deltaE_for_substitution
	(
			int alternate_state,
			core::PackerEnergy & prev_node_energy
	);
	void commit_considered_substitution();

	// <directed_design>
	/// @brief returns the change in weighted energy that would be induced
	/// by switching this node from its current state into another state
	void project_deltaE_for_substitution
	(
			int alternate_state,
			core::PackerEnergy & deltaE_unweighted,
			core::PackerEnergy & prevE_unweighted,
			core::PackerEnergy& deltaE_weighted,
			core::PackerEnergy& prevE_weighted,
			ObjexxFCL::FArray2D< core::PackerEnergy > const& weights
	);
	core::PackerEnergy get_weighted_energy_with_higher_indexed_nodes(ObjexxFCL::FArray2D< core::PackerEnergy > const& weights) const;
	// </directed_design>

	inline
	void acknowledge_neighbors_state_substitution(
			int edge_to_altered_neighbor,
			core::PackerEnergy new_edge_energy,
			int other_node_new_state,
			SparseMatrixIndex const & other_node_new_state_sparse_info);

	SparseMatrixIndex const &
	get_sparse_mat_info_for_state(int state) const;

	SparseMatrixIndex const &
	get_sparse_mat_info_for_curr_state() const;


	int  get_num_states_for_aa_type(int aa_type) const;
	void print_internal_energies() const;

	void update_internal_energy_sums();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;


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
	void get_absent_rots( ObjexxFCL::FArray1_int & rots_absent );

	int get_num_states_in_file();
	int & get_aatypes_for_file_states();

	int & get_aatypes_for_states();
	int & get_num_file_states_for_aa();
	int & get_file_states_2_instance_states_array();

	bool get_node_corresponded_to_file_node();
	*/

protected:
	void update_internal_vectors();

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

	void
	set_alt_aa_offsets_from_edge(int edge_index, ObjexxFCL::FArray2D_int const & offsets);

	inline
	PDEdge const * get_incident_pd_edge( int index ) const;

	inline
	PDEdge * get_incident_pd_edge( int index );

	inline
	PDInteractionGraph const * get_pdig_owner() const;

	inline
	PDInteractionGraph * get_pdig_owner();

	int num_aa_types_;
	utility::vector1< int > num_states_for_aatype_;
	std::vector< SparseMatrixIndex > sparse_mat_info_for_state_;
	std::vector< core::PackerEnergy > one_body_energies_;

	ObjexxFCL::FArray3D_int aa_offsets_for_edges_;
	ObjexxFCL::FArray2D_int num_states_for_aa_type_for_higher_indexed_neighbor_;
	std::vector< int > neighbors_curr_state_;
	std::vector< SparseMatrixIndex > neighbors_curr_state_sparse_info_;
	std::vector< ObjexxFCL::FArray1A< core::PackerEnergy > > edge_matrix_ptrs_;


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

/*
	//variables for I/O
	int num_states_in_file_;
	ObjexxFCL::FArray1D_int instance_states_2_file_states_;
	ObjexxFCL::FArray1D_int file_states_2_instance_states_;
	ObjexxFCL::FArray1D_int aa_types_for_file_states_;
	ObjexxFCL::FArray1D_int aa_types_for_instance_states_; //useful only for I/O
	ObjexxFCL::FArray1D_int num_file_states_for_aa_;
*/
	bool alternate_state_is_being_considered_;

	//no default constructor, uncopyable
	PDNode();
	PDNode( PDNode const & );
	PDNode & operator = ( PDNode const & );
};

class PDEdge : public PrecomputedPairEnergiesEdge
{
public:
	PDEdge(InteractionGraphBase* owner, int first_node_ind, int second_node_ind);
	virtual ~PDEdge();
	virtual void set_sparse_aa_info(ObjexxFCL::FArray2_bool const & sparse_conn_info);
	virtual void force_aa_neighbors(int node1aa, int node2aa);
	virtual void force_all_aa_neighbors();
	virtual bool get_sparse_aa_info( int node1aa, int node2aa) const;
	virtual void add_to_two_body_energy(int const, int const, core::PackerEnergy const);
	virtual void
	add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array );
	virtual
	void set_two_body_energy(int const, int const, core::PackerEnergy const);
	virtual
	void clear_two_body_energy(int const, int const);
	virtual core::PackerEnergy get_two_body_energy( int const, int const ) const;

	virtual void declare_energies_final();
	virtual void prepare_for_simulated_annealing();
	//virtual unsigned int getMemoryUsageInBytes() const;

	core::PackerEnergy get_current_two_body_energy();

	void acknowledge_state_change(
			int node_ind,
			int new_state,
			SparseMatrixIndex const & new_state_sparse_info,
			core::PackerEnergy & new_energy
	);

	/// @brief updates bookkeeping information when one of the two nodes enters its
	/// "unassigned" state.
	void acknowledge_state_zeroed( int node_ind );

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

	inline void acknowledge_substitution(
			int substituted_node_index,
			core::PackerEnergy const curr_state_energy,
			int nodes_new_state,
			SparseMatrixIndex const & nodes_new_state_sparse_info
	);

	/// @brief Returns the array of offsets into the sparse two-body energy table
	/// for amino-acid neighbors.  Used in transferring information from edges
	/// onto nodes for cache efficiency.
	ObjexxFCL::FArray2D_int const & get_offsets_for_aatypes( );

	/// @brief returns an FArray of the number of states for each amino acid type for the
	/// higher-indexed node
	utility::vector1< int > const & get_second_node_num_states_per_aa();

	int get_two_body_table_size() const;
	core::PackerEnergy & get_edge_table_ptr();

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	ObjexxFCL::FArray2D< core::PackerEnergy >
	get_aa_submatrix_energies(
		int node1aa,
		int node2aa
	) const;

/*
	void read_edge_energies_from_file( std::ifstream & infile );
	static void skip_over_edge_energies_from_file
	(
			std::ifstream & infile,
			int num_aa,
			ObjexxFCL::FArray1_int & num_file_states_for_aa_node1,
			ObjexxFCL::FArray1_int & num_file_states_for_aa_node2
	);
	void write_edge_energies_to_file( std::ofstream & outfile );
*/

	virtual void set_edge_weight( Real weight );

protected:

	//Hooks for SASAEdge< V, E, G > class
	void declare_energies_final_no_deletion();
	void prepare_for_simulated_annealing_no_deletion();
	bool pd_edge_table_all_zeros() const;

private:

	inline
	PDNode const * get_pd_node( int index ) const;

	inline
	PDNode * get_pd_node( int index );

	inline
	PDInteractionGraph const * get_pdig_owner() const;

	inline
	PDInteractionGraph * get_pdig_owner();

	void drop_small_submatrices_where_possible( core::PackerEnergy epsilon );
	void drop_zero_submatrices_where_possible();

private: // Data

	AminoAcidNeighborSparseMatrix< core::PackerEnergy > two_body_energies_;
	core::PackerEnergy curr_state_energy_;
	bool energies_updated_since_last_prep_for_simA_;

	//no default constructor, uncopyable
	PDEdge();
	PDEdge( PDEdge const & );
	PDEdge & operator = ( PDEdge const & );
};

class PDInteractionGraph : public PrecomputedPairEnergiesInteractionGraph
{
public:
	PDInteractionGraph(int num_nodes);
	virtual void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

	virtual core::PackerEnergy get_one_body_energy_for_node_state( int node, int state);

	/// @brief sets the number of amino acid types present.
	//virtual void set_num_aatypes(int);
	virtual int  get_num_aatypes() const;

	virtual void add_edge(int node1, int node2);

	virtual void blanket_assign_state_0();
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state);
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states);
	virtual void consider_substitution
	(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	);
	virtual core::PackerEnergy commit_considered_substitution();
	/// @ brief O(1) total energy report.  Protected read access for derived classes.
	virtual core::PackerEnergy get_energy_current_state_assignment();

	/// @brief older scheme for memory accounting -- replace this asap
	/// @brief returns the number of floats used in all edge two-body energy tables
	virtual int get_edge_memory_usage() const;
	/// @brief outputs the current state for each node, useful for debugging
	virtual void print_current_state_assignment() const;
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy ) {};

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;


/*
	//Methods for I/O
	void prepare_to_read_energies_from_file();
	void declare_finished_reading_from_file();
	void set_num_file_aatypes( int num_file_aatypes );
	int get_num_file_aatypes();
	void set_num_nodes_in_file( int num_nodes_in_file );
	void set_node_correspondence( int instance_node, int file_node );
	void set_num_states_for_file_node(int node, int num_file_states);
	void set_aa_for_file_node_state( int file_node, int file_state, int state_aa );
	void set_correspondence_for_state(int node, int state, int file_state);
	int get_correspondence_for_state(int node, int state );
	bool get_node_corresponded_to_file_node( int node );

	int get_num_rots_absent_from_file(int node);
	void get_absent_rots(int node, ObjexxFCL::FArray1_int & rots_absent );

	void read_edge_energies_from_file( std::ifstream & infile );
	void write_edge_energies_to_file( std::ofstream & outfile );
*/

	/// @brief a user may define subsets of the vertex set for which they would like to
	/// know the internal energy sum.
	virtual core::PackerEnergy get_energy_sum_for_vertex_group( int group_id );

	// <directed_design>
	core::PackerEnergy get_weighted_energy(ObjexxFCL::FArray2D< core::PackerEnergy > const &weights) const;
	core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states, ObjexxFCL::FArray2D< core::PackerEnergy > const& weights);
	virtual void consider_substitution
	(
	 int node_ind,
	 int new_state,
	 core::PackerEnergy & deltaE_unweighted,
	 core::PackerEnergy & prevE_unweighted, // !!! remember this is the energy just for this node
	 core::PackerEnergy & deltaE_weighted,
	 core::PackerEnergy & prevE_weighted,
	 ObjexxFCL::FArray2D< core::PackerEnergy > const& weights
	 );
	virtual core::PackerEnergy commit_considered_substitution(ObjexxFCL::FArray2D< core::PackerEnergy > const& weights);
	// </directed_design>

	/// @brief Override the InteractionGraphBase class's implementation of this function
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

protected:
	virtual unsigned int getMemoryUsageInBytes() const;

	/// @brief factory method that instantiates a PDNode.
	virtual NodeBase* create_new_node( int node_index, int num_states);
	/// @brief factory method that instantiates a PDEdge
	virtual EdgeBase* create_new_edge( int index1, int index2);

	//Hooks for SASAInterationGraph< V, E, G >
	core::PackerEnergy get_energy_PD_current_state_assignment();
	/// @brief removes numerical drift that can accumulate over the course of
	/// many state assignment changes within simulated annealing
	void update_internal_energy_totals();

	inline
	PDNode* get_pd_node(int index) const
	{       return (PDNode*) get_node( index );}

	inline
	PDEdge const * get_pd_edge( int node1, int node2 ) const
	{
		return (PDEdge const *) find_edge( node1, node2 );
	}

	inline
	PDEdge * get_pd_edge( int node1, int node2 )
	{
		return (PDEdge *) find_edge( node1, node2 );
	}


private:
	int num_aa_types_;
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;


	//variables for I/O
	int num_nodes_in_file_;
	int num_file_aa_types_;
	ObjexxFCL::FArray1D_int file_node_2_instance_node_;
	ObjexxFCL::FArray1D_int instance_node_2_file_node_;
	ObjexxFCL::FArray1D< ObjexxFCL::FArray1D_int > aa_types_for_states_on_file_nodes_;
	ObjexxFCL::FArray1D< ObjexxFCL::FArray1D_int > num_file_states_for_aa_for_node_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	PDInteractionGraph();
	PDInteractionGraph( PDInteractionGraph const & );
	PDInteractionGraph & operator = ( PDInteractionGraph const & );
};

inline
PDEdge const *
PDNode::get_incident_pd_edge( int index ) const
{
	return static_cast< PDEdge const * > ( get_incident_edge( index ));
}

inline
PDEdge *
PDNode::get_incident_pd_edge( int index )
{
	return static_cast< PDEdge * > ( get_incident_edge( index ));
}

inline
PDInteractionGraph const *
PDNode::get_pdig_owner() const
{
	return static_cast< PDInteractionGraph const * > (get_owner());
}

inline
PDInteractionGraph *
PDNode::get_pdig_owner()
{
	return static_cast< PDInteractionGraph * > (get_owner());
}


inline
PDNode const *
PDEdge::get_pd_node( int index ) const
{
	return static_cast< PDNode const *> (get_node( index ));
}

inline
PDNode *
PDEdge::get_pd_node( int index )
{
	return static_cast< PDNode *> (get_node( index ));
}

inline
PDInteractionGraph const *
PDEdge::get_pdig_owner() const
{
	return static_cast< PDInteractionGraph const * > (get_owner());
}

inline
PDInteractionGraph *
PDEdge::get_pdig_owner()
{
	return static_cast< PDInteractionGraph * > (get_owner());
}

/// @brief static method that looks up the two body energy when the
/// node with the smaller index on an edge is considering an alternate state
///
/// @param first_node_alt_state - [in] - the alternate state for the lower-indexed node
/// @param second_node_orig_state - [in] - the current state for the higher-indexed node
/// @param second_node_orig_state_sparse_info - [in] - the sparse matrix info for
///    the higher-indexed node
/// @param first_node_state_offset_minus_1 - [in] - part of the sparse matrix info
///   for the lower-indexed node where 1 is subtracted from the state offset.
/// @param second_node_num_states_per_aatype - [in] - number of states with current aa
///   type for node 2
/// @param aa_neighbor_offset - [in] - offset for the amino-acid neighbor pair for
///    the sparse two-body energy table
/// @param edge_energy_table - [in] - the proxy FArray pointing at the edge table
///   connecting the two nodes.
///
inline
float
PDEdge::get_alternate_state_energy_first_node(
	int first_node_alt_state,
	int second_node_orig_state,
	SparseMatrixIndex const & second_node_orig_state_sparse_info,
	int first_node_state_offset_minus_1,
	int second_node_curr_num_states_per_aatype,
	int aa_neighbor_offset,
	ObjexxFCL::FArray1< core::PackerEnergy > & edge_energy_table
)
{

	if (first_node_alt_state == 0 || second_node_orig_state == 0) {
		return 0.0f;
	} else {
		return AminoAcidNeighborSparseMatrix< float >::get(
			second_node_orig_state_sparse_info,
			first_node_state_offset_minus_1,
			second_node_curr_num_states_per_aatype,
			aa_neighbor_offset,
			edge_energy_table );
	}
}

/// @brief update bookkeeping information when one of the nodes an edge is incident
/// upon changes state
///
/// @param substituted_node_index - [in] - index of the node that chagned its state
/// @param curr_state_energy - [in] - the two body energy given the new state
/// @param nodes_new_state - [in] - the state the node just transitioned into
/// @param nodes_new_state_sparse_info - [in] - sparse matrix info for the new state
///
inline
void
PDEdge::acknowledge_substitution(
	int substituted_node_index,
	float const curr_state_energy,
	int nodes_new_state,
	SparseMatrixIndex const & nodes_new_state_sparse_info
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = curr_state_energy;

	get_pd_node( node_not_substituted )->
	acknowledge_neighbors_state_substitution (
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state,
		nodes_new_state_sparse_info
	);

	return;
}

/// @brief static method that looks up the two body energy when the
/// node with the larger index on an edge is considering an alternate state
///
/// @param first_node_orig_state - [in] - the current state for the lower-indexed node
/// @param second_node_alt_state - [in] - the alt state for the higher-indexed node
/// @param first_node_orig_state_sparse_info - [in] - the sparse matrix info for
///   the lower-indexed node
/// @param second_node_alt_state_sparse_info - [in] - the sparse matrix info for
///   the higher-indexed node
/// @param first_node_state_offset_minus_1 - [in] - part of the sparse matrix info
///   for the lower-indexed node where 1 is subtracted from the state offset.
/// @param second_node_alt_state_num_states_per_aatype - [in] - number of states
///   with alternate aa type for node 2
/// @param aa_neighbor_offset - [in] - offset for the amino-acid neighbor pair for
///   the sparse two-body energy table
/// @param edge_energy_table - [in] - the proxy FArray pointing at the edge table
///   connecting the two nodes.
///
inline
float
PDEdge::get_alternate_state_energy_second_node(
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
		return 0.0f;
	} else {
		return AminoAcidNeighborSparseMatrix< float >::get(
			first_node_orig_state_sparse_info,
			second_node_alternate_state_sparse_info,
			second_node_alt_state_num_states_per_aatype,
			aa_neighbor_offset,
			edge_energy_table );
	}
}

/// @brief updates bookkeeping arrays for when a neighbor has changed its state
///
/// @param edge_to_altered_neighbor - [in] - the index for the edge that connects
/// 	this node to the node that just changed its state
/// @param new_edge_energ - [in] - the pair energy between this node in its current
///	state and the new state of the node that just changed its state
/// @param other_node_new_state - [in] - the state the neighbor just adopted
/// @param other_node_new_state_sparse_info - [in] - the sparse-matrix info
///	corresponding to the neighbor's new state
///
inline
void PDNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	float new_edge_energy,
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

/// @detais iterates across the incident edges for a node in two phases:
/// in the first phase, it examines edges leading to higher-indexed nodes
/// in the second phase, it examines edges leading to smaller-indexed nodes.
/// for cache efficiency, all of the amino-acid-neighbor-offset information
/// that each edge calculates is stored on the nodes themselves.  The edges
/// are never touched; rather, their private information is stored on the nodes
/// and handed to static member functions of the PDEdge class.  This "store
/// edge information on the nodes" strategy gives me performance equivalent
/// to the previous energy2b lookup tables.
///
/// @param alternate_state - [in] - the alternate state to consider
/// @param previous_energy_for_node - [out] - the old energy1b/energy2b sum for this
/// node; used by simulate annealing.
inline
float
PDNode::project_deltaE_for_substitution(
	int alternate_state,
	float & prev_energy_for_node
)
{

	alternate_state_is_being_considered_ = true;
	//std::cout << "proj_deltaE: node -  " << get_node_index() << " alt state " << alternate_state << "...";

	alternate_state_ = alternate_state;
	alt_state_sparse_mat_info_ = sparse_mat_info_for_state_[ alternate_state];
	alternate_state_one_body_energy_ = one_body_energies_[ alternate_state ];
	//std::cout << "alternate_state_one_body_energy_: " << alternate_state_one_body_energy_ << std::endl;
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
	prev_energy_for_node = curr_state_total_energy_;

	int alt_state_num_states_per_aa_type =
		num_states_for_aatype_[ alt_state_sparse_mat_info_.get_aa_type() ];
	int alt_state_for_aa_type_minus_1 =
		alt_state_sparse_mat_info_.get_state_ind_for_this_aa_type() - 1;
	int nstates_offset =
		num_states_for_aa_type_for_higher_indexed_neighbor_.index(1,1) - 1;
	int aa_neighb_linear_index_offset = aa_offsets_for_edges_.
		index(1, 1, alt_state_sparse_mat_info_.get_aa_type() ) - 1;


	for (int ii = 1; ii <= get_num_edges_to_smaller_indexed_nodes();
		++ii, aa_neighb_linear_index_offset += num_aa_types_) {

		alternate_state_two_body_energies_[ ii ] =
			get_incident_pd_edge(ii)->
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
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		//if ( alternate_state_two_body_energies_[ ii ] != 0.0 ) {
		//	std::cout << "( " << get_index_of_adjacent_node( ii ) << " , " << alternate_state_two_body_energies_[ ii ] <<" ) ";
		//}
	}

	for (int ii = get_num_edges_to_smaller_indexed_nodes() + 1;
		ii <= get_num_incident_edges();
		++ii, aa_neighb_linear_index_offset += num_aa_types_,
		nstates_offset += num_aa_types_) {

		alternate_state_two_body_energies_[ ii ] =
			get_incident_pd_edge(ii)->
				get_alternate_state_energy_first_node(
				alternate_state_,
				neighbors_curr_state_[ii],
				//alt_state_sparse_mat_info_,
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
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		//if ( alternate_state_two_body_energies_[ ii ] != 0.0 ) {
		//	std::cout << "( " << get_index_of_adjacent_node( ii ) << " , " << alternate_state_two_body_energies_[ ii ] <<" ) ";
		//}
	}

	//std::cerr<< "..done" << std::endl;

	return alternate_state_total_energy_ - curr_state_total_energy_;

}


} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif //SPARSE_PD_INTERACTION_GRAPH_H
