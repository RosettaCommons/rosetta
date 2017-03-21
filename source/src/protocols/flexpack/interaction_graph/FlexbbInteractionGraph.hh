// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/FlexbbIteractionGraph.hh
/// @brief  Declaration for flexible-backbone-packing interaction graph interface & base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_FlexbbInteractionGraph_hh
#define INCLUDED_protocols_flexpack_interaction_graph_FlexbbInteractionGraph_hh

/// Unit headers
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.fwd.hh>

/// Package headers
#include <protocols/flexpack/interaction_graph/FlexbbSparseMatrixIndex.hh>

/// Project headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

/// Utility headers

// ObjexxFCL Headers

/// C++ headers
/// TEMP

#include <utility/vector0_bool.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

class FlexbbNode : public core::pack::interaction_graph::NodeBase
{
public:
	typedef core::pack::interaction_graph::NodeBase parent;
	typedef core::PackerEnergy PackerEnergy;
	typedef core::Size Size;

public:

	FlexbbNode( FlexbbInteractionGraph * owner, int node_id, int num_states);
	virtual ~FlexbbNode();
	virtual void print() const;

	void set_num_distinct_backbones( int nbbconfs );
	int  get_num_distinct_backbones() const { return num_bb_; }
	void set_num_states_per_backbone( utility::vector1< int > const & );
	int  get_bb_for_state( int  state ) const { return state_info_[ state ].get_bb(); }
	int  get_num_states_for_bb( int  bbconf ) const;
	int  get_state_offset_for_bb( int  bbconf ) const;
	void get_states_on_curr_bb(
		utility::vector1< Size > & state_list,
		int offset
	) const;
	void get_all_states(
		utility::vector1< Size > & rotlist,
		int offset
	) const;

	void set_closest_states_on_other_bbs( ObjexxFCL::FArray2D_int const & );

	void set_amino_acid_types( utility::vector1< int > const & );
	ObjexxFCL::FArray1A_int getNumStatesPerAAPerBB( int aa );


	virtual void add_to_one_body_energies( ObjexxFCL::FArray1< PackerEnergy > & energies );
	virtual void add_to_one_body_energy( int state, PackerEnergy energy );
	virtual void update_one_body_energy( int state, PackerEnergy energy);
	virtual void zero_one_body_energies();
	virtual bool state_unassigned() const;

	PackerEnergy get_one_body_energy( int state ) const;

	virtual void prepare_for_simulated_annealing();

	int get_current_state() const;
	int get_backbone_for_current_state() const;
	void write_current_state_to_state_array( ObjexxFCL::FArray1_int & nodes_states);

	/// @brief Preliminatry bookkeeping for the node the Graph contacted about a backbone move.
	/// This node is the "root" of the DFS traversal.
	//void
	//register_contacted_node_for_bb_jump();

	/// @brief Preliminary bookkeeping for DFS-backbone-motion inducing state substitution
	/// initiated by a separate node.  Returns false if the current state is unassigned.
	//bool prepare_for_bb_jump( int alt_bb );
	bool inform_edges_of_alt_state_before_bbjump();

	/// @brief Did the alterenate state assigned during a flexbb move actually preserve the currently assigned bb?
	//bool
	//bb_move_actually_kept_original_bb() const;

	bool
	state_has_same_backbone_as_current( int state ) const
	{
		return state_info_[ state ].get_bb() == current_state_info_.get_bb();
	}


	//inline PackerEnergy project_deltaE_for_alt_bb_state (
	// int alternate_state,
	// PackerEnergy & curr_frag_etotal,
	// bool & valid_motion
	//);

	//inline PackerEnergy project_deltaE_for_alt_bb(
	// int alternate_backbone,
	// PackerEnergy & curr_frag_etotal,
	// bool & valid_motion
	//);

	//PackerEnergy get_altE_for_bb_move
	//(
	// PackerEnergy & curr_frag_etotal
	//);

	//void commit_considered_substitution();
	//void commit_considered_substitution( FArray1_int & rotamer_on_node );
	//void commit_alt_bb_substitution( FArray1_int & rotamer_on_node );

	//inline
	//void acknowledge_neighbors_state_substitution(
	// int edge_to_altered_neighbor,
	// PackerEnergy new_edge_energy
	//);

	//virtual void note_last_considered_bbsub_uncommitted() = 0;

	int  get_num_states_for_aa_type_for_bb( int aa_type, int bb ) const;
	void print_internal_energies();

	void update_internal_energy_sums();

	virtual unsigned int count_dynamic_memory() const;

protected:

	/// Downcast pointers to incident edges, adjacent nodes, and the owning graph
	///                                                    ^--- Oxford comma. See Vampire Weekend's objection.
	inline
	FlexbbEdge const * get_incident_flexbb_edge( int index ) const;

	inline
	FlexbbEdge * get_incident_flexbb_edge( int index );

	inline
	FlexbbNode const * get_adjacent_flexbb_node( int index ) const;

	inline
	FlexbbNode * get_adjacent_flexbb_node( int index );

	inline
	FlexbbInteractionGraph const * get_flexbbig_owner() const;

	inline
	FlexbbInteractionGraph * get_flexbbig_owner();

protected:
	void update_internal_vectors();
	void inform_edges_considered_fixedbb_substition_uncommitted();

public:
	/// Read access to private data to both derived classes and the world
	int num_aa_types() const { return num_aa_types_; }
	int num_bbconfs() const { return num_bb_; }
	utility::vector1< int > const & num_states_for_bb() { return num_states_for_bb_; }
	utility::vector1< int > const & state_offsets_for_bb() { return state_offsets_for_bb_; }
	ObjexxFCL::FArray2D_int const & num_states_for_aa_for_bb() { return num_states_for_aa_for_bb_; }
	ObjexxFCL::FArray2D_int const & state_offsets_for_aa_for_bb() { return state_offsets_for_aa_for_bb_; }

	ObjexxFCL::FArray2D_int const & closest_state_on_alt_bb() const { return closest_state_on_alt_bb_; }
	FlexbbSparseMatrixIndex const & state_info( int state ) const { return state_info_[ state ]; }
	utility::vector1< PackerEnergy > const & one_body_energies() const { return one_body_energies_; }

	int current_state() const { return current_state_; }
	FlexbbSparseMatrixIndex const & current_state_info() const { return current_state_info_; }
	int curr_bb() const { return current_state_info_.get_bb(); }
	PackerEnergy curr_state_one_body_energy() const { return curr_state_one_body_energy_; }
	PackerEnergy curr_state_total_energy() const { return curr_state_total_energy_; }
	utility::vector1< PackerEnergy > const & curr_state_two_body_energies() { return curr_state_two_body_energies_; }


	int alternate_state() const { return alternate_state_; }
	FlexbbSparseMatrixIndex const & alternate_state_info() const { return alternate_state_info_; }
	int alt_bb() const { return alternate_state_info_.get_bb(); }


	PackerEnergy alternate_state_one_body_energy() const { return alternate_state_one_body_energy_;}
	PackerEnergy alternate_state_total_energy() const { return alternate_state_total_energy_; }
	utility::vector1< PackerEnergy > const & alternate_state_two_body_energies() const { return alternate_state_two_body_energies_; }
	PackerEnergy alternate_state_two_body_energies( int ind ) const { return alternate_state_two_body_energies_[ ind ]; }

	utility::vector1< bool > const & edge_connects_flexsegmate() const { return edge_connects_flexsegmate_; }

	/// @brief Convention for acumulating energies between flexsegmates during a backbone move:
	/// The smaller-indexed node counts the interaction, the larger-indexed node does not.
	bool count_energy_to_node_in_my_fragtotalE( int edge_index ) const { return edge_index > get_num_edges_to_smaller_indexed_nodes(); }
	bool alternate_state_is_being_considered() const { return alternate_state_is_being_considered_; }

protected:

	/// Limited write acces to private data
	void set_considering_alternate_state() { alternate_state_is_being_considered_ = true; }
	void set_current_state( int setting ) { current_state_ = setting; current_state_info_ = state_info_[ setting ]; }
	void set_curr_state_one_body_energy( PackerEnergy setting ) { curr_state_one_body_energy_ = setting; }
	void set_curr_state_total_energy( PackerEnergy setting ) { curr_state_total_energy_ = setting; }
	void inc_curr_state_total_energy( PackerEnergy setting ) { curr_state_total_energy_ += setting; }
	void set_curr_state_two_body_energies( Size index, PackerEnergy setting ) { curr_state_two_body_energies_[ index ] = setting; }

	void set_alternate_state( int setting ) { alternate_state_ = setting; alternate_state_info_ = state_info_[ setting ]; }
	void set_alternate_state_one_body_energy( PackerEnergy setting ) {
		//if ( get_node_index() == 17 && alternate_state_ == 63 ) { std::cout << "SETTING ALT STATE 1BE: " << setting << std::endl; }
		alternate_state_one_body_energy_ = setting;
	}
	void set_alternate_state_total_energy( PackerEnergy setting ) { alternate_state_total_energy_ = setting; }
	void inc_alternate_state_total_energy( PackerEnergy setting ) { alternate_state_total_energy_ += setting; }
	void set_alternate_state_two_body_energies( Size index, PackerEnergy setting ) { alternate_state_two_body_energies_[ index ] = setting;; }

	void reset_all_rotamer_substitution_bookkeeping_data();

	void partial_state_assignment( int new_state );
	void inform_incident_edges_about_partial_state_assignment();

	//// @brief bookkeeping for DFS traversal.  Derived classes should consult this at the beginning of node
	/// "get_energies" traversals.  Calling this function for a node effectively declares that the code which
	/// follows will project the delta energy: after the first time this functio is called, it will return true.
	/// Must be proceeded by a call to prepare_for_bb_jump().
	/// Depricated.
	///bool energies_already_projected() const;

	/// @brief Copy over the energy data local to this node from "alternate_*" to "current_*"
	void copy_alternate_to_current();

	/// @brief Have all incident edges copy their data from "alternate_*" to "current_*"
	void have_edges_copy_alternate_to_current();

	/// @brief Have a subset of the incident edges copy their data from "alternate_*" to "current_*"
	/// where this subset includes all edges to non-flexsegmates as well as all
	/// all upper edges to flexsegmates.  Threadsafe if it can be fruitfully parallelized.
	void have_edges_copy_alternate_to_current_following_flexbb_accept();

private:

	int num_aa_types_;
	int num_bb_;

	utility::vector1< int > num_states_for_bb_;
	utility::vector1< int > state_offsets_for_bb_;
	ObjexxFCL::FArray2D_int num_states_for_aa_for_bb_; // Indexed (aa, bb)
	ObjexxFCL::FArray2D_int state_offsets_for_aa_for_bb_; // Indexed (aa, bb)
	ObjexxFCL::FArray2D_int closest_state_on_alt_bb_;

	/// vector0 so that assigning state 0 does not index out-of-bounds.
	utility::vector0< FlexbbSparseMatrixIndex > state_info_; // stores aa index and bb index for each state.
	//utility::vector1< int > num_states_for_aatype_; // ?!

	utility::vector1< PackerEnergy > one_body_energies_;

	int current_state_;
	FlexbbSparseMatrixIndex current_state_info_;
	PackerEnergy curr_state_one_body_energy_;
	PackerEnergy curr_state_total_energy_;
	utility::vector1< PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
	FlexbbSparseMatrixIndex alternate_state_info_;
	PackerEnergy alternate_state_one_body_energy_;
	PackerEnergy alternate_state_total_energy_;
	utility::vector1< PackerEnergy > alternate_state_two_body_energies_;

	utility::vector1< bool > edge_connects_flexsegmate_;

	bool alternate_state_is_being_considered_;
	//bool told_edges_alt_state_for_bb_move_;
	//mutable bool projected_energies_for_bb_move_; // is it necessary that the calling function is const?
	//bool resolved_considered_bb_move_;
	//bool node_contacted_by_graph_about_bb_move_;

	// uncopyable
	FlexbbNode();
	FlexbbNode( FlexbbNode const & );

};

class FlexbbEdge : public core::pack::interaction_graph::EdgeBase
{
public:
	typedef core::pack::interaction_graph::EdgeBase parent;
	typedef core::PackerEnergy PackerEnergy;
	typedef core::Real Real;
	typedef core::Size Size;

public:
	FlexbbEdge(
		FlexbbInteractionGraph * owner,
		int first_node_ind,
		int second_node_ind
	);

	virtual ~FlexbbEdge();

	//virtual void set_nodes_from_same_flexseg( bool same_flexseg );
	bool get_nodes_from_same_flexseg() const {
		return nodes_part_of_same_flexseg_;
	}

	/// @brief Called by FlexbbNode in prepare_for_bb_jump: Edges must know the alternate
	/// states that nodes are considering.  Precondition: alt_state_ for a fixed node must
	/// match its cur_state_
	void
	set_alt_state( int node_index, int new_state, FlexbbSparseMatrixIndex const & state_info );

	void
	acknowledge_partial_state_assignment( int node_index, int new_state, FlexbbSparseMatrixIndex const & state_info );

	/// @brief Copy alt data to current data after a state substitution
	void
	note_state_substitution_accepted();

	/// @brief After a rejected state substitution, the FlexbbNode will call this function
	/// to reset the alt_state data to establish the invariant that the alternate state held
	/// on the flexbb edges reflect the the current state of those nodes not considering a
	/// state substitution
	void
	reset_alternate_states_for_uncommited_substitution();

	virtual unsigned int count_dynamic_memory() const;

	PackerEnergy cur_energy() const { return cur_energy_; }
	PackerEnergy alt_energy() const { debug_assert( alt_e_up_to_date_ ); return alt_energy_; }


protected:

	inline
	FlexbbNode const * get_flexbb_node( int index ) const
	{ return static_cast< FlexbbNode const * > (get_node( index )); }

	inline
	FlexbbNode * get_flexbb_node( int index )
	{ return static_cast< FlexbbNode * > (get_node( index )); }

	inline
	FlexbbInteractionGraph const * get_flexbbig_owner() const;

	inline
	FlexbbInteractionGraph * get_flexbbig_owner();

protected:

	inline
	int
	num_bb( int node ) {
		debug_assert( node == 0 || node == 1 );
		return nodes_num_bb_[ node ];
	}

	inline
	bool
	nodes_part_of_same_flexseg() const {
		return nodes_part_of_same_flexseg_;
	}

	inline
	bool
	nodes_considering_bb_move() const {
		return nodes_considering_bb_move_;
	}

	inline
	void
	set_nodes_considering_bb_move( bool status ) {
		nodes_considering_bb_move_ = status;
	}

	int nodes_cur_state( int node ) const { debug_assert( node == 0 || node == 1 ); return nodes_cur_state_[ node ]; }
	FlexbbSparseMatrixIndex const & nodes_cur_info( int node ) const    { debug_assert( node == 0 || node == 1 ); return  nodes_cur_info_[ node ]; }
	int nodes_alt_state( int node ) const { debug_assert( node == 0 || node == 1 ); return  nodes_alt_state_[ node ]; }
	FlexbbSparseMatrixIndex const & nodes_alt_info( int node ) const    { debug_assert( node == 0 || node == 1 ); return  nodes_alt_info_[ node ]; }
	bool alt_e_up_to_date() const { return alt_e_up_to_date_; }


	/// @brief Set the currently assigned state for a node; node == 0 || 1
	void set_nodes_cur_state( int node, int setting ) { nodes_cur_state_[ node ] = setting; }
	void set_nodes_cur_info( int node, FlexbbSparseMatrixIndex const & setting ) { nodes_cur_info_[ node ] = setting; }
	/// @brief Set the altnernate state being considered for a node; node == 0 || 1
	void set_nodes_alt_state( int node, int setting ) { nodes_alt_state_[ node ] = setting; }
	void set_nodes_alt_info( int node, FlexbbSparseMatrixIndex const & setting ) { nodes_alt_info_[ node ] = setting; }
	void set_cur_energy( PackerEnergy setting ) { cur_energy_ = setting; }
	void set_alt_energy( PackerEnergy setting ) { alt_energy_ = setting; alt_e_up_to_date_ = true; }

	void copy_alternate_to_current();
	void set_node_state_to_zero( int which_node );

private:

	bool const nodes_part_of_same_flexseg_;
	int  nodes_num_bb_[2];

	int                     nodes_cur_state_[ 2 ];
	FlexbbSparseMatrixIndex nodes_cur_info_[ 2 ];

	int                     nodes_alt_state_[ 2 ];
	FlexbbSparseMatrixIndex nodes_alt_info_[ 2 ];

	PackerEnergy cur_energy_;
	PackerEnergy alt_energy_;

	bool alt_e_up_to_date_;
	bool nodes_considering_bb_move_;


	// uncopyable
	FlexbbEdge();
	FlexbbEdge( FlexbbEdge const & );
};


class FlexbbInteractionGraph : public core::pack::interaction_graph::InteractionGraphBase
{
public:
	typedef core::pack::interaction_graph::InteractionGraphBase parent;
	typedef core::PackerEnergy PackerEnergy;
	typedef core::pack::interaction_graph::EdgeBase EdgeBase;

public:
	enum Subsitution {SC_ONLY, BOTH_SC_AND_BB};

public:
	virtual ~FlexbbInteractionGraph();
	FlexbbInteractionGraph(int num_nodes);

	virtual void initialize(core::pack::rotamer_set::RotamerSetsBase const & rot_sets );

private:
	/// Private functions called during initialize()
	void set_num_flexsegs(int num_flexsegs);
	void set_total_num_backbones( int num_backbones );

public:
	int  get_num_aa_types() const { debug_assert( num_aa_types_ != 0 ); return num_aa_types_; }
	//void set_representitive_node_for_flexseg( int flexseg, int node_index);
	//int  get_flexseg_representative( int flexseg ) const { return flexseg_representative_[ flexseg ]; }
	utility::vector1< Size > const &
	flexseg_members( int flexseg ) const {
		return flexseg_members_[ flexseg ];
	}
	int  get_flexseg_for_bb( int bb ) const { return flexseg_for_bb_[ bb ]; }
	int  get_flexseg_bb_offset( int flexseg_id ) const { return flexseg_bb_offset_[ flexseg_id ]; }
	bool nodes_from_same_flexseg( int node1, int node2 ) const;

private:
	void set_num_bb_for_node( int node, int numbb);
	void set_num_states_per_backbone_for_node( int node, utility::vector1< int > const & states_per_bb);

public:
	int  get_num_states_per_backbone_for_node( int node, int bb ) const;
	int  get_bb_for_state( int node, int state ) const;

private:
	void set_aatypes_for_node(int node_ind, utility::vector1< int > const & aatypes);
	void set_closest_states_on_other_bbs( int node_index, ObjexxFCL::FArray2D_int const & );

public:
	//void set_edge_connecting_nodes_on_same_flexseg(int, int); // knowable at edge construction time!

	virtual
	void consider_backbone_move(
		int bb_id,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_flexseg_energy,
		bool & valid_motion,
		int & num_nodes_changing_state
	) = 0;

	virtual
	void consider_bbmove_w_state_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_flexseg,
		bool & valid_motion,
		int & num_nodes_changing_state
	) = 0;

	virtual
	PackerEnergy
	commit_considered_backbone_move(
		ObjexxFCL::FArray1_int & rotamer_on_node
	) = 0;

	void get_accessible_states(
		Subsitution move_mode,
		utility::vector1< Size > & rotlist
	) const;

	void get_backbone_list(
		utility::vector1< Size > & bblist
	) const;

	/// @brief Is the backbone conformation (in the global enumertion of backbone conformations) already
	/// assigned to the network?  False if any residue on the flexible segment that this bbid corresponds to
	/// is assigned state 0.
	bool get_backbone_currently_assigned( int bbid ) const;

	/// @brief FlexbbNodes will ask: am I allowed to have a state that breaks the backbone?
	/// There are brief periods when the backbone is "broken" as the graph assigns new states to
	/// nodes on the same flexible segment.
	bool get_enforce_bb_contiguity() const;

	/// @brief Owner keeps a count of the number of nodes undergoing a
	/// simultaneous rotamer substitution as the backbone moves.
	void increment_count_nodes_in_flexseg();

	virtual unsigned int count_dynamic_memory() const;

protected:

	/// Downcasts
	FlexbbNode const * get_flexbb_node( int index ) const
	{ return static_cast< FlexbbNode const * > (get_node( index )); }

	FlexbbNode * get_flexbb_node( int index )
	{ return static_cast< FlexbbNode * > (get_node( index )); }

	FlexbbEdge const * find_flexbb_edge( int node1, int node2 ) const
	{
		core::pack::interaction_graph::EdgeBase const * edge = find_edge( node1, node2 );
		if ( edge ) return static_cast< FlexbbEdge const * > ( edge );
		else return 0;
	}

	FlexbbEdge * find_flexbb_edge( int node1, int node2 )
	{
		core::pack::interaction_graph::EdgeBase * edge = find_edge( node1, node2 );
		if ( edge ) return static_cast< FlexbbEdge * > ( edge );
		else return 0;
	}

	FlexbbEdge const * cast_flexbb_edge( EdgeBase const * edge ) const
	{ debug_assert( mine( edge ) ); return static_cast< FlexbbEdge const * > ( edge ); }

	FlexbbEdge * cast_flexbb_edge( EdgeBase * edge )
	{ debug_assert( mine( edge ) ); return static_cast< FlexbbEdge * > ( edge ); }


protected:

	void set_enforce_bb_contiguity(bool);
	void note_bbjump_substitution() { last_sub_attempted_backbone_move_ = true; }
	void note_fixedbb_substitution() { last_sub_attempted_backbone_move_ = false; last_considered_fixedbb_sub_unresolved_ = true; }
	bool last_considered_substitution_kept_backbone_fixed() const { return !last_sub_attempted_backbone_move_; }
	bool last_considered_substitution_moved_the_backbone() const { return last_sub_attempted_backbone_move_; }

	bool last_considered_backbone_sub_unresolved() const { return last_considered_backbone_sub_unresolved_; }
	bool last_considered_substitution_unresolved() const {
		return last_considered_backbone_sub_unresolved_ || last_considered_fixedbb_sub_unresolved_;
	}
	void note_last_considered_substitution_resolved();

	/// @details at the beginning of backbone-changing substitutions,
	/// derived classes must invoke this function to get a proper count
	/// of the number of nodes undergoing a simultaneous rotamer substitution
	void reset_node_in_moving_flexseg_count();
	int get_num_nodes_changing_state() const;

	PackerEnergy total_energy_current_state_assignment() const { return total_energy_current_state_assignment_; }
	PackerEnergy total_energy_alternate_state_assignment() const { return total_energy_alternate_state_assignment_; }

	/// @brief Only allowed to ask for the node considering an alternate
	/// state during a fixed-backbone substitution
	int node_considering_alt_state() const { debug_assert( node_considering_alt_state_ != 0 ); return node_considering_alt_state_; }

	/// @brief Only allowed to ask for the flexible segment considering an alternate
	/// backbone conformation during a backbone-moving substitution
	int flexseg_considering_alt_bb() const { debug_assert( flexseg_considering_alt_bb_ != 0 ); return flexseg_considering_alt_bb_; }

	bool last_considered_backbone_sub_valid() const { return last_considered_backbone_sub_valid_; }

	int flexseg_for_moltenres( int moltenres ) const {
		return flexseg_for_moltenres_[ moltenres ];
	}

	void set_total_energy_current_state_assignment( PackerEnergy setting );
	void set_total_energy_alternate_state_assignment( PackerEnergy setting ) { total_energy_alternate_state_assignment_ = setting; }

	/// @brief Track the last node at which a fixed-backbone substitution took place
	void set_node_considering_alt_state( int setting ) {
		flexseg_considering_alt_bb_ = 0;
		node_considering_alt_state_ = setting;
	}
	/// @brief Track the last flexible segment at which a moving-backbone substitution took place
	void set_flexseg_considering_alt_bb( int setting ) {
		flexseg_considering_alt_bb_ = setting;
		node_considering_alt_state_ = 0;
	}
	void set_last_considered_backbone_sub_valid( bool setting ) { last_considered_backbone_sub_valid_ = setting; }


	virtual void update_internal_energy_totals();

private:

	int num_aa_types_;

	PackerEnergy total_energy_current_state_assignment_;
	PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;
	int flexseg_considering_alt_bb_;

	int num_flexible_segments_;
	int num_total_bb_;
	//utility::vector1< int > flexseg_representative_;
	utility::vector1< utility::vector1< Size > > flexseg_members_;
	utility::vector1< int > num_bb_alternatives_for_flexseg_;
	utility::vector1< int > flexseg_for_bb_;
	utility::vector1< int > flexseg_bb_offset_;
	utility::vector1< int > flexseg_for_moltenres_;

	bool enforce_bb_contiguity_;
	bool last_sub_attempted_backbone_move_;
	bool last_considered_backbone_sub_valid_;
	bool last_considered_backbone_sub_unresolved_;
	bool last_considered_fixedbb_sub_unresolved_;
	int  num_nodes_changing_state_;


	int num_commits_since_last_update_;
	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 20; //1024; // 2^10

	FlexbbInteractionGraph();
	FlexbbInteractionGraph( FlexbbInteractionGraph const & );


};

inline
FlexbbEdge const * FlexbbNode::get_incident_flexbb_edge( int index ) const
{ return static_cast< FlexbbEdge const * > ( get_incident_edge( index )); }

inline
FlexbbEdge * FlexbbNode::get_incident_flexbb_edge( int index )
{ return static_cast< FlexbbEdge * > ( get_incident_edge( index )); }

inline
FlexbbNode const * FlexbbNode::get_adjacent_flexbb_node( int index ) const
{ return static_cast< FlexbbNode const * > ( get_adjacent_node( index )); }

inline
FlexbbNode * FlexbbNode::get_adjacent_flexbb_node( int index )
{ return static_cast< FlexbbNode * > ( get_adjacent_node( index )); }


inline
FlexbbInteractionGraph const * FlexbbNode::get_flexbbig_owner() const
{ return static_cast< FlexbbInteractionGraph const * > (get_owner()); }

inline
FlexbbInteractionGraph * FlexbbNode::get_flexbbig_owner()
{ return static_cast< FlexbbInteractionGraph * > (get_owner()); }


inline
FlexbbInteractionGraph const * FlexbbEdge::get_flexbbig_owner() const
{ return static_cast< FlexbbInteractionGraph const * > (get_owner()); }

inline
FlexbbInteractionGraph * FlexbbEdge::get_flexbbig_owner()
{ return static_cast< FlexbbInteractionGraph * > (get_owner()); }


}
}
}

#endif

