// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/InteractionGraphBase.hh
/// @brief  Interaction graph base class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_InteractionGraphBase_hh
#define INCLUDED_core_pack_interaction_graph_InteractionGraphBase_hh

// Unit Headers
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>


// Package Headers

#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// STL Headers
#include <iosfwd>
// AUTO-REMOVED #include <vector>
#include <list>
// AUTO-REMOVED #include <assert.h>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1_bool.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class NodeBase
{
public:
	virtual ~NodeBase();
	NodeBase( InteractionGraphBase*, int node_id, int num_states);
	int get_num_states() const;
	std::list< EdgeBase* >::iterator add_edge(EdgeBase* edge_ptr);
	void drop_edge(std::list< EdgeBase* >::iterator edge_iterator);
	void drop_all_edges();
	EdgeBase* find_edge(int other_node_index) const;

	virtual void assign_zero_state() = 0;
	virtual void prepare_for_simulated_annealing() = 0;
	virtual void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energy1b ) = 0;
	virtual void add_to_one_body_energy( int state, core::PackerEnergy energy ) = 0;
	virtual void update_one_body_energy( int state, core::PackerEnergy energy) = 0;
	virtual void zero_one_body_energies() = 0;
	virtual void print() const = 0;
	virtual bool state_unassigned() const = 0;

	void depth_first_connected_component_counting();

	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

	std::list< EdgeBase * >::const_iterator edge_list_begin();
	std::list< EdgeBase * >::const_iterator edge_list_end();

protected:
	void update_edge_vector();

public:

	//Read access to private data granted to derived classes
	//inlined for efficiency
	inline
	int get_node_index() const
	{
		return node_index_;
	}

	inline
	int get_num_incident_edges() const
	{
		return num_incident_edges_;
	}

	inline
	int get_num_edges_to_smaller_indexed_nodes() const
	{
		return num_edges_to_smaller_indexed_nodes_;
	}

	inline
	int get_num_edges_to_larger_indexed_nodes() const
	{
		return num_edges_to_larger_indexed_nodes_;
	}

public:

	/// These functions are public for the sake of writing good unit tests

	inline
	EdgeBase const * get_incident_edge( int index ) const
	{
		assert( edge_vector_up_to_date_ );
		return incident_edge_vector_[ index ];
	}

	inline
	EdgeBase * get_incident_edge( int index )
	{
		assert( edge_vector_up_to_date_ );
		return incident_edge_vector_[ index ];
	}

	inline
	int get_index_of_adjacent_node( int index ) const
	{
		assert( edge_vector_up_to_date_ );
		return adjacent_node_ind_[ index ];
	}

	inline
	NodeBase const * get_adjacent_node( int index ) const
	{
		assert( edge_vector_up_to_date_ );
		return adjacent_node_[ index ];
	}


	inline
	NodeBase * get_adjacent_node( int index )
	{
		assert( edge_vector_up_to_date_ );
		return adjacent_node_[ index ];
	}

protected:
	inline
	bool get_edge_vector_up_to_date() const
	{
		return edge_vector_up_to_date_;
	}

	inline
	InteractionGraphBase const * get_owner() const
	{
		return owner_;
	}

	inline
	InteractionGraphBase * get_owner()
	{
		return owner_;
	}


private:

	int node_index_;
	int num_states_;
	int num_incident_edges_;
	int num_edges_to_smaller_indexed_nodes_;
	int num_edges_to_larger_indexed_nodes_;
	std::list< EdgeBase* > incident_edge_list_;
	std::vector< EdgeBase* > incident_edge_vector_;
	std::vector< int > adjacent_node_ind_;
	std::vector< NodeBase* > adjacent_node_;
	bool edge_vector_up_to_date_;
	InteractionGraphBase* owner_;

	//no default constructor, uncopyable
	NodeBase();
	NodeBase( NodeBase const & );
	NodeBase & operator = ( NodeBase & );
};

class EdgeBase
{
public:
	virtual ~EdgeBase();
	EdgeBase(InteractionGraphBase* owner, int first_node_ind, int second_node_ind);

	int get_other_ind(int callers_index) const;
	NodeBase* get_other_node(int callers_index) const;
	int get_first_node_ind() const;
	int get_second_node_ind() const;
	void set_pos_in_owners_list( std::list< EdgeBase* >::iterator edge_iterator);
	void set_pos_in_node_edgevector(int callers_index, int position);

	bool same_edge(int node1, int node2) const;

	virtual void declare_energies_final() = 0;
	virtual void prepare_for_simulated_annealing() = 0;

	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

	Real edge_weight() const {
		return edge_weight_;
	}

	virtual void set_edge_weight( Real weight ) = 0;

protected:

	//Read access to private data granted to derived classes
	inline
	int get_node_index( int index ) const
	{
		assert( index == 0 || index == 1 );
		return node_indices_[ index ];
	}

	inline
	int get_num_states_for_node( int index ) const
	{
		assert( index == 0 || index == 1 );
		return num_node_states_[ index ];
	}

	inline
	NodeBase const *
	get_node( int index ) const
	{
		assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}

	inline
	NodeBase *
	get_node( int index )
	{
		assert( index == 0 || index == 1 );
		return nodes_[ index ];
	}

	inline
	int get_edges_position_in_nodes_edge_vector( int index ) const
	{
		assert( index == 0 || index == 1 );
		return pos_in_nodes_edge_vector_[ index ];
	}

public:
	inline
	InteractionGraphBase const * get_owner() const
	{
		return owner_;
	}

	inline
	InteractionGraphBase * get_owner()
	{
		return owner_;
	}

protected:

	/// @brief is a node the first or second node this edge is incident upon?
	inline
	int
	which_node( int node_index ) const {
		assert( node_index == node_indices_[ 0 ] || node_index == node_indices_[ 1 ] );
		return ( node_index == node_indices_[ 0 ] ? 0 : 1 );
	}

	/// @brief protected setter of the edge weight.  To be called by derived
	/// classes after they have completed the conversion from the previous edge weighting
	/// to the new edge weighting.
	void
	edge_weight( Real );

private:
	int node_indices_[2];
	int num_node_states_[2];
	NodeBase* nodes_[2];
	int pos_in_nodes_edge_vector_[2];
	std::list< EdgeBase* >::iterator pos_in_nodes_edge_list_[2];
	std::list< EdgeBase* >::iterator pos_in_owners_edge_list_;
	InteractionGraphBase* owner_;

	/// Allow the arbitrary scaling of energies for each edge.
	/// The derived classes have the responsibility of scaling each
	/// energy as it is added to the edge, and also updating all of
	/// the edge energies if the edge weight changes *after* all
	/// the edge energies have been stored -- that is, by dividing each
	/// energy by the old weight and multiplying by the new weight.
	Real edge_weight_;


	//no default constructor, uncopyable
	EdgeBase();
	EdgeBase( EdgeBase const & );
	EdgeBase & operator = ( EdgeBase & );

};

class InteractionGraphBase : public utility::pointer::ReferenceCount
{
public:
	virtual ~InteractionGraphBase();

	InteractionGraphBase(int num_nodes);

	virtual int get_num_nodes_v() const
	{
		return get_num_nodes();
	}

	inline
	int get_num_nodes() const
	{
		return num_ig_nodes_;
	}

	virtual void initialize( rotamer_set::RotamerSetsBase const & rot_sets ) = 0;

	void set_num_states_for_node(int node, int num_states);
	int  get_num_states_for_node(int node) const;
	int  get_num_total_states() const {return num_total_states_;}
	virtual core::PackerEnergy get_one_body_energy_for_node_state( int node, int state) = 0;
	void add_edge( int node1, int node2);
	bool get_edge_exists(int node1, int node2);
	void drop_all_edges_for_node( int node );

	void print_vertices() const;
	virtual void print() const {};
	void output_connectivity(std::ostream & os) const;
	void output_dimacs(std::ostream & os) const;

	/// @brief iterate across edges and nodes and allow them to prepare
	/// for simulated annealing
	virtual void prepare_for_simulated_annealing();

	virtual void  blanket_assign_state_0() = 0;
	virtual core::PackerEnergy set_state_for_node(int node_ind, int new_state) = 0;
	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states) = 0;
	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node) = 0;
	virtual core::PackerEnergy commit_considered_substitution() = 0;
	virtual core::PackerEnergy get_energy_current_state_assignment() = 0;

	void set_edge_weight( int node1, int node2, Real edge_weight );
	Real get_edge_weight( int node1, int node2 ) const;

	virtual int get_edge_memory_usage() const = 0;
	virtual void print_current_state_assignment() const = 0;
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE ) = 0;

	bool any_vertex_state_unassigned() const;

	void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & one_body_energies );
	void update_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & old_energy1b, ObjexxFCL::FArray1< core::PackerEnergy > & new_energy1b);
	void zero_one_body_energies_for_node( int node );

	void add_to_nodes_one_body_energy
	(
		int node_ind,
		utility::vector1< core::PackerEnergy > const & one_body_energies
	);

	void add_to_nodes_one_body_energy
	(
		int node_ind,
		ObjexxFCL::FArray1< core::PackerEnergy > const & one_body_energies
	);

	/// @brief interface to PrecomputedPairEnergiesNode::add_to_nodes_one_body_energy
	void add_to_nodes_one_body_energy
	(
		int node_ind,
		int state_id,
		core::PackerEnergy const one_body_energy
	);


	void set_number_of_energy_sum_vertex_groups( int num_groups );
	void set_vertex_member_of_group( int vertex, int group );
	void print_vertex_groups();
	virtual core::PackerEnergy get_energy_sum_for_vertex_group( int group_id ) = 0;
	int count_connected_components_and_initialize_vertex_groups();
	void note_vertex_reached( int node_index );
	bool vertex_already_reached( int node_index );

	inline
	bool get_vertex_member_of_energy_sum_group( int node_index, int group_id )
	{
		assert(  num_energy_sum_groups_ != -1 &&
			node_index > 0 && node_index <= num_ig_nodes_ &&
			group_id > 0 && group_id <= num_energy_sum_groups_);
		return energy_sum_group_membership_( node_index, group_id );
	}

	virtual unsigned int getTotalMemoryUsage() const;

	/// Methods for iterating across edges of the interaction graph.
	/// Protected access since the raw edge lists contain non-const pointers.

	/// @brief set the Graph's (single) edge list iterator to the beginning of the edge list
	/// for a particular node
	void reset_edge_list_iterator_for_node( int node_index ) const;
	/// @brief increment the (single) edge list iterator to the next element
	void increment_edge_list_iterator() const;
	/// @brief test: have we arrived at the edge list end?
	bool edge_list_iterator_at_end() const;
	/// @brief return a const reference to an edge pointed at by the list iterator
	EdgeBase const & get_edge() const;

	friend class NodeBase;
	friend class EdgeBase;

protected:

	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

	void drop_edge(std::list< EdgeBase* >::iterator edge);

public:

	/// The following functions provide access to the nodes and edges in the graph
	/// though, their use is strongly discouraged except for in writing unit tests
	/// to ensure that the graphs are properly implemented.

	EdgeBase const * find_edge(int node1, int node2) const;
	EdgeBase * find_edge(int node1, int node2);

	virtual NodeBase* create_new_node( int node_index, int num_states) = 0;
	virtual EdgeBase* create_new_edge( int index1, int index2) = 0;

	inline
	NodeBase* get_node( int index ) const
	{
		assert( index > 0 && index <= num_ig_nodes_ );
		return ig_nodes_[ index ];
	}

	inline
	int get_num_edges() const
	{
		return ig_edge_list_.size();
	}

	inline
	std::list< EdgeBase* >::iterator get_edge_list_begin()
	{
		return ig_edge_list_.begin();
	}

	inline
	std::list< EdgeBase* >::const_iterator get_edge_list_begin() const
	{
		return ig_edge_list_.begin();
	}

	inline
	std::list< EdgeBase* >::const_iterator get_edge_list_end() const
	{
		return ig_edge_list_.end();
	}

protected:

	inline
	int get_node_state_offset( int index ) const
	{
		assert( index > 0 && index <= num_ig_nodes_ );
		return node_state_offsets_[ index ];
	}

	bool
	mine( EdgeBase const * edge ) const {
		return edge->get_owner() == this;
	}

private:
	int num_ig_nodes_;
	std::vector< NodeBase* > ig_nodes_;
	std::list< EdgeBase* > ig_edge_list_;

	std::vector< int > node_state_offsets_;
	int num_total_states_;
	mutable EdgeBase* focused_edge_;
	mutable std::list< EdgeBase * >::const_iterator focused_edge_iterator_;
	mutable std::list< EdgeBase * >::const_iterator focused_edge_iterator_end_;

	int num_energy_sum_groups_;
	ObjexxFCL::FArray2D_bool energy_sum_group_membership_;
	ObjexxFCL::FArray1D_int component_membership_;

	//no default constructor, uncopyable
	InteractionGraphBase();
	InteractionGraphBase( InteractionGraphBase const &);
	InteractionGraphBase & operator = (InteractionGraphBase const & );
};

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif //INTERACTION_GRAPH_BASE_CLASS_H
