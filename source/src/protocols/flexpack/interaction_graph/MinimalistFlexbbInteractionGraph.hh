// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/MinimalistFlexbbInteractionGraph.hh
/// @brief  Class declaration for minimimalist on-the-fly RPE calculating FlexbbInteractionGraph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_MinimalistFlexbbInteractionGraph_hh
#define INCLUDED_protocols_flexpack_interaction_graph_MinimalistFlexbbInteractionGraph_hh

/// Unit headers
#include <protocols/flexpack/interaction_graph/MinimalistFlexbbInteractionGraph.fwd.hh>

/// Package headers
#include <protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

class MinimalistFlexbbNode : public OTFFlexbbNode
{
public:
	typedef OTFFlexbbNode parent;
	typedef FlexbbNode    grandparent;

public:

	MinimalistFlexbbNode( MinimalistFlexbbInteractionGraph * owner, int node_id, int num_states);
	~MinimalistFlexbbNode() override;

	/// Virtual functions from NodeBase not covered by FlexbbNode
	void assign_zero_state() override;
	void prepare_for_simulated_annealing() override;
	void print() const override;

	/// Virtual functions from OTFFlexbbNode

	/// Other functions

	PackerEnergy
	project_deltaE_for_substitution(
		int alternate_state,
		PackerEnergy & prev_energy_for_node
	);

	/// @brief Traverse incident edges and inform them of the alternate
	/// state being considered at this node.  Return true if the alternate
	/// state is valid.
	bool
	prepare_for_altbb_move_to_state( int alt_state );

	/// @brief Set the alternate state to be the one on the alternate backbone
	/// most similar to the currently assigned state.
	/// Traverse incident edges and inform them of the alternate
	/// state being considered at this node.  Return true if the alternate
	/// state is valid.
	bool
	prepare_for_altbb_move_to_closest_state( int alt_bb );

	PackerEnergy
	get_frag_energy_for_alt_bb_state();

	PackerEnergy
	get_frag_energy_for_curr_bb_state_and_finalize_alt_energy_total();

	/*
	PackerEnergy
	project_deltaE_with_backbone_move(
	int new_state,
	PackerEnergy & prev_energy_for_flexseg,
	bool & valid_motion
	);

	PackerEnergy
	project_deltaE_for_backbone_move(
	int alt_bb,
	PackerEnergy & prev_energy_for_flexseg,
	bool & valid_motion
	);
	*/

	void
	commit_considered_substitution();

	void
	commit_considered_substitution( ObjexxFCL::FArray1_int & state_on_node );

	void
	commit_alt_bb_substitution( ObjexxFCL::FArray1_int & state_on_node );

	void
	acknowledge_neighbors_substitution(
		int which_edge,
		PackerEnergy alternate_energy
	);

	void
	resolve_uncommitted_substitution();

	/// @brief assign a new state to this node and return the change in energy induced
	PackerEnergy
	assign_state( int new_state );

	/// @brief For assinging network states to the interaction graph: when substituting
	/// two adjacent nodes i and j from states r and s to states r' and s', this
	/// function, along with its partner function, complete_partial_state_assignment,
	/// guarantees that neither the energies r' with s nor s' with r are computed
	/// which would be wasteful.  This function is meant to be called by the
	/// owning interaction graph and not meant to be called by the outside world.
	void
	partially_assign_state( int new_state );

	/// @brief See description of partially_assign_state.
	void
	complete_partial_state_assignment();

	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

protected:

	/// Downcast pointers to incident edges, adjacent nodes, and the owning graph
	inline
	MinimalistFlexbbEdge const * get_incident_minimalistflexbb_edge( int index ) const;

	inline
	MinimalistFlexbbEdge * get_incident_minimalistflexbb_edge( int index );

	inline
	MinimalistFlexbbNode const * get_adjacent_minimalistflexbb_node( int index ) const;

	inline
	MinimalistFlexbbNode * get_adjacent_minimalistflexbb_node( int index );

	inline
	MinimalistFlexbbInteractionGraph const * get_minimalistflexbbig_owner() const;

	inline
	MinimalistFlexbbInteractionGraph * get_minimalistflexbbig_owner();


protected:

	PackerEnergy
	get_altE_for_bb_move( PackerEnergy & curr_frag_etotal );


private:


};

class MinimalistFlexbbEdge : public OTFFlexbbEdge
{
public:
	typedef OTFFlexbbEdge parent;
	typedef core::Real Real;

public:
	MinimalistFlexbbEdge( MinimalistFlexbbInteractionGraph * owner, int node1, int node2 );
	~MinimalistFlexbbEdge() override;

	/// Virtual functions from EdgeBase
	void declare_energies_final() override;
	void prepare_for_simulated_annealing() override;
	void set_edge_weight( Real weight ) override;


	PackerEnergy
	get_alt_stateE();

	void
	acknowledge_substitution( int node_changing_state );

	void
	acknowledge_state_zeroed( int node_index );

	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

protected:

	/// Downcasts

	inline
	MinimalistFlexbbNode const * get_minimalistflexbb_node( int index ) const;

	inline
	MinimalistFlexbbNode * get_minimalistflexbb_node( int index );

	inline
	MinimalistFlexbbInteractionGraph const * get_minimalistflexbbig_owner() const;

	inline
	MinimalistFlexbbInteractionGraph * get_minimalistflexbbig_owner();

private:

};

class MinimalistFlexbbInteractionGraph : public OTFFlexbbInteractionGraph
{
public:
	typedef OTFFlexbbInteractionGraph parent;

public:
	MinimalistFlexbbInteractionGraph( int num_nodes );
	~MinimalistFlexbbInteractionGraph() override;

	/// Virtual functions from InteractionGraphBase
	void initialize( core::pack_basic::RotamerSetsBase const & rot_sets ) override;
	PackerEnergy get_one_body_energy_for_node_state( int node, int state) override;
	void  blanket_assign_state_0() override;
	PackerEnergy set_state_for_node(int node_ind, int new_state) override;
	PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states) override;
	void consider_substitution(
		int node_ind,
		int new_state,
		PackerEnergy & delta_energy,
		PackerEnergy & prev_energy_for_node) override;
	PackerEnergy commit_considered_substitution() override;
	PackerEnergy get_energy_current_state_assignment() override;
	int get_edge_memory_usage() const override;
	void print_current_state_assignment() const override;
	void set_errorfull_deltaE_threshold( PackerEnergy deltaE ) override;
	PackerEnergy get_energy_sum_for_vertex_group( int group_id ) override;


	/// Virtual functions from FlexbbInteractionGraph
	void consider_backbone_move(
		int bb_id,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_flexseg_energy,
		bool & valid_motion,
		int & num_nodes_changing_state
	) override;

	void consider_bbmove_w_state_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_flexseg,
		bool & valid_motion,
		int & num_nodes_changing_state
	) override;

protected:
	void
	complete_deltaE_prediction_for_bbmove(
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_flexseg,
		bool & valid_motion,
		int & num_nodes_changing_state
	);

private:

	PackerEnergy
	commit_considered_backbone_move(
		ObjexxFCL::FArray1_int & rotamer_on_node
	) override;

	/// Virtual functions from OTFFlexbbInteractionGraph


	/// @brief returns the change in energy that would be induced by switching this node
	/// from its current state into another state on the same backbone
	PackerEnergy
	project_deltaE_for_substitution(
		int alternate_state,
		PackerEnergy & prev_energy_for_node
	);

protected:
	/// Downcasts
	inline
	MinimalistFlexbbNode const * get_minimalistflexbb_node( int index ) const;

	inline
	MinimalistFlexbbNode * get_minimalistflexbb_node( int index );

	inline
	MinimalistFlexbbEdge const * find_minimalist_flexbb_edge( int node1, int node2 ) const;

	inline
	MinimalistFlexbbEdge * find_minimalist_flexbb_edge( int node1, int node2 );

	inline
	MinimalistFlexbbEdge const * cast_minimalist_flexbb_edge( EdgeBase const * edge ) const;

	inline
	MinimalistFlexbbEdge * cast_minimalist_flexbb_edge( EdgeBase * edge );

protected:
	/// Virtual functions from InteractionGraphBase
	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

	core::pack::interaction_graph::NodeBase * create_new_node( int node_index, int num_states) override;
	core::pack::interaction_graph::EdgeBase * create_new_edge( int index1, int index2) override;


	void resolve_uncommitted_substitution();

	/// Virtual functions from OTFFlexbbInteractionGraph

private:


};

/// Node Downcast
inline
MinimalistFlexbbEdge const *
MinimalistFlexbbNode::get_incident_minimalistflexbb_edge( int index ) const
{ return static_cast< MinimalistFlexbbEdge const * > ( get_incident_edge( index )); }

inline
MinimalistFlexbbEdge *
MinimalistFlexbbNode::get_incident_minimalistflexbb_edge( int index )
{ return static_cast< MinimalistFlexbbEdge * > ( get_incident_edge( index )); }

inline
MinimalistFlexbbNode const *
MinimalistFlexbbNode::get_adjacent_minimalistflexbb_node( int index ) const
{ return static_cast< MinimalistFlexbbNode const * > ( get_adjacent_node( index )); }

inline
MinimalistFlexbbNode *
MinimalistFlexbbNode::get_adjacent_minimalistflexbb_node( int index )
{ return static_cast< MinimalistFlexbbNode * > ( get_adjacent_node( index )); }


inline
MinimalistFlexbbInteractionGraph const *
MinimalistFlexbbNode::get_minimalistflexbbig_owner() const
{ return static_cast< MinimalistFlexbbInteractionGraph const * > (get_owner()); }

inline
MinimalistFlexbbInteractionGraph *
MinimalistFlexbbNode::get_minimalistflexbbig_owner()
{ return static_cast< MinimalistFlexbbInteractionGraph * > (get_owner()); }

/// Edge Downcasts

inline
MinimalistFlexbbNode const *
MinimalistFlexbbEdge::get_minimalistflexbb_node( int index ) const
{ return static_cast< MinimalistFlexbbNode const * > (get_node( index )); }

inline
MinimalistFlexbbNode *
MinimalistFlexbbEdge::get_minimalistflexbb_node( int index )
{ return static_cast< MinimalistFlexbbNode * > (get_node( index )); }

inline
MinimalistFlexbbInteractionGraph const *
MinimalistFlexbbEdge::get_minimalistflexbbig_owner() const
{ return static_cast< MinimalistFlexbbInteractionGraph const * > (get_owner()); }

inline
MinimalistFlexbbInteractionGraph *
MinimalistFlexbbEdge::get_minimalistflexbbig_owner()
{ return static_cast< MinimalistFlexbbInteractionGraph * > (get_owner()); }


/// Graph Downcasts
inline
MinimalistFlexbbNode const *
MinimalistFlexbbInteractionGraph::get_minimalistflexbb_node( int index ) const
{ return static_cast< MinimalistFlexbbNode const * > (get_node( index )); }

inline
MinimalistFlexbbNode *
MinimalistFlexbbInteractionGraph::get_minimalistflexbb_node( int index )
{ return static_cast< MinimalistFlexbbNode * > (get_node( index )); }

inline
MinimalistFlexbbEdge const *
MinimalistFlexbbInteractionGraph::find_minimalist_flexbb_edge( int node1, int node2 ) const
{
	core::pack::interaction_graph::EdgeBase const * edge = find_edge( node1, node2 );
	if ( edge ) return static_cast< MinimalistFlexbbEdge const * > ( edge );
	else return nullptr;
}

inline
MinimalistFlexbbEdge *
MinimalistFlexbbInteractionGraph::find_minimalist_flexbb_edge( int node1, int node2 )
{
	core::pack::interaction_graph::EdgeBase * edge = find_edge( node1, node2 );
	if ( edge ) return static_cast< MinimalistFlexbbEdge * > ( edge );
	else return nullptr;
}

inline
MinimalistFlexbbEdge const *
MinimalistFlexbbInteractionGraph::cast_minimalist_flexbb_edge( EdgeBase const * edge ) const
{ debug_assert( mine( edge ) ); return static_cast< MinimalistFlexbbEdge const * > ( edge ); }

inline
MinimalistFlexbbEdge *
MinimalistFlexbbInteractionGraph::cast_minimalist_flexbb_edge( EdgeBase * edge )
{ debug_assert( mine( edge ) ); return static_cast< MinimalistFlexbbEdge * > ( edge ); }


}
}
}

#endif

