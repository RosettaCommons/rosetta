// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/FixedBBInteractionGraph.hh
/// @brief  Interaction graph base class for fixed-backbone packing
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_FixedBBInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_FixedBBInteractionGraph_hh

// Unit Headers
#include <core/pack/interaction_graph/FixedBBInteractionGraph.fwd.hh>

// Package Headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

/// This layer of abstraction between the InteractionGraphBase and the
/// various fixed-backbone interaction graphs is primarily to allow
/// outside users to READ from an interaction graph -- to grab data
/// that the interaction graph stores.  This is not so much for the
/// writing of information -- in particular, edge energies -- to an
/// interaction graph.

namespace core {
namespace pack {
namespace interaction_graph {

class FixedBBNode : public NodeBase
{
public:
	virtual ~FixedBBNode() {}

	FixedBBNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states)
	:
		NodeBase( owner, node_id, num_states )
	{}

};

class FixedBBEdge : public EdgeBase
{
public:
	virtual ~FixedBBEdge() {}

	FixedBBEdge(
		InteractionGraphBase* owner,
		int first_node_ind,
		int second_node_ind)
	:
		EdgeBase( owner, first_node_ind, second_node_ind )
	{}

	virtual
	void set_sparse_aa_info(ObjexxFCL::FArray2_bool const & sparse_conn_info) = 0;

	virtual
	bool get_sparse_aa_info( int node1aa, int node2aa) const = 0;

	virtual
	void force_aa_neighbors(int node1aa, int node2aa) = 0;

	virtual
	void force_all_aa_neighbors() = 0;

	virtual core::PackerEnergy get_two_body_energy( int const, int const ) const = 0;

};

class FixedBBInteractionGraph : public InteractionGraphBase
{
public:
	virtual ~FixedBBInteractionGraph() {}
	FixedBBInteractionGraph( int num_nodes )
	:
		InteractionGraphBase( num_nodes )
	{}

	//virtual void set_num_aatypes( int ) = 0;
	virtual int get_num_aatypes() const = 0;

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

	/// @brief interface to PDEdge::set_sparse_aa_info
	void set_sparse_aa_info_for_edge
	(
		int node1,
		int node2,
		ObjexxFCL::FArray2_bool const & sparse_conn_info
	);

	/// @brief returns true if node1aa and node2aa are amino acid neighbors
	bool get_sparse_aa_info_for_edge
	(
		int node1,
		int node2,
		int node1aa,
		int node2aa
	);

	/// @brief interface to FixedBBEdge::force_aa_neighbors
	void force_aa_neighbors_for_edge
	(
		int node1,
		int node2,
		int node1aa,
		int node2aa
	);

	/// @brief interface to PDEdge::force_aa_neighbors
	void force_all_aa_neighbors_for_edge
	(
		int node1,
		int node2
	);


	/// @brief interface to FixedBBEdge::get_two_body_energy
	///  - returns the state pair energy
	virtual core::PackerEnergy get_two_body_energy_for_edge
	(
		int node1,
		int node2,
		int state_node1,
		int state_node2
	) const;

protected:
	/// Downcasts

	inline FixedBBNode const * get_fixedbb_node( int node_index ) const;
	inline FixedBBNode       * get_fixedbb_node( int node_index );
	inline FixedBBEdge const * get_fixedbb_edge( int node1, int node2 ) const;
	inline FixedBBEdge       * get_fixedbb_edge( int node1, int node2 );

};

inline
FixedBBNode const * FixedBBInteractionGraph::get_fixedbb_node( int node_index ) const
{
	return static_cast< FixedBBNode const * > (get_node( node_index ));
}

inline
FixedBBNode       * FixedBBInteractionGraph::get_fixedbb_node( int node_index )
{
	return static_cast< FixedBBNode * > (get_node( node_index ));
}

inline
FixedBBEdge const * FixedBBInteractionGraph::get_fixedbb_edge( int node1, int node2 ) const
{
	return static_cast< FixedBBEdge const * > (find_edge( node1, node2 ));
}

inline
FixedBBEdge       * FixedBBInteractionGraph::get_fixedbb_edge( int node1, int node2 )
{
	return static_cast< FixedBBEdge * > (find_edge( node1, node2 ));
}

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif
