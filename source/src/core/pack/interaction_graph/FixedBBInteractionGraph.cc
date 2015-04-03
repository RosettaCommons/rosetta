// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FixedBBInteractionGraph.cc
/// @brief  Precomputed interaction graph class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Rosetta Headers
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/exit.hh>

// STL Headers
#include <iostream>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {

/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param sparse_conn_info - [in] - the boolean 2-D array of amino-acid neighbor info
void FixedBBInteractionGraph::set_sparse_aa_info_for_edge(
	int node1,
	int node2,
	FArray2_bool const & sparse_conn_info)
{
	FixedBBEdge* edge = get_fixedbb_edge( node1, node2 );
	if (edge == NULL)
	{	std::cerr <<
			"WARNING:: you've input sparse aa info for an edge that does not exist"
			<< std::endl;
		return;
	}
	edge->set_sparse_aa_info( sparse_conn_info );
	return;
}


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param aa_node1 - [in] - the amino acid type for the node with the smaller index
/// @param aa_node2 - [in] - the amino acid type for the node with the larger index
bool FixedBBInteractionGraph::get_sparse_aa_info_for_edge(
	int node1,
	int node2,
	int node1aa,
	int node2aa
)
{
	FixedBBEdge* edge = get_fixedbb_edge( node1, node2 );
	if (edge == NULL) {
		std::cerr << "WARNING:: you've requested sparse aa info for an edge that does not exist" << std::endl;
		return false;
	}

	return edge->get_sparse_aa_info( node1aa, node2aa );
}


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param aa_node1 - [in] - the amino acid type for the node with the smaller index
/// @param aa_node2 - [in] - the amino acid type for the node with the larger index
void FixedBBInteractionGraph::force_aa_neighbors_for_edge
(
	int node1,
	int node2,
	int aa_node1,
	int aa_node2
)
{
	FixedBBEdge* edge = get_fixedbb_edge( node1, node2 );
	if (edge == NULL) {
		return;
	}
	edge->force_aa_neighbors(aa_node1, aa_node2);
	return;
}


/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param aa_node1 - [in] - the amino acid type for the node with the smaller index
/// @param aa_node2 - [in] - the amino acid type for the node with the larger index
void FixedBBInteractionGraph::force_all_aa_neighbors_for_edge
(
	int node1,
	int node2
)
{
	FixedBBEdge* edge = get_fixedbb_edge( node1, node2 );
	if (edge == NULL) {
		add_edge( node1, node2 );
		edge = (FixedBBEdge*) find_edge(node1, node2);

		FArray2D_bool all_aa_neighbors( get_num_aatypes() , get_num_aatypes(), true);
		edge->set_sparse_aa_info( all_aa_neighbors );
	} else {
		edge->force_all_aa_neighbors();
	}
	return;
}

/// @param node1 - [in] - the index of the smaller-indexed node
/// @param node2 - [in] - the index of the larger-indexed node
/// @param state_node1 - [in] - state on smaller-indexed node
/// @param state_node2 - [in] - state on larger-indexed node
core::PackerEnergy
FixedBBInteractionGraph::get_two_body_energy_for_edge(
	int node1,
	int node2,
	int state_node1,
	int state_node2
) const
{
	FixedBBEdge const * edge = get_fixedbb_edge( node1, node2 );
	if (edge == NULL) {
		return 0;
	}
	return edge->get_two_body_energy( state_node1, state_node2 );
}

/// @details Base class provides no implementation for this functionality
bool
FixedBBInteractionGraph::aa_submatrix_energies_retrievable() const
{
	return false;
}

/// @details Derived classes wishing to respond to the get_aa_submatrix_energies_for_edge
/// must implement this function as well, which allows the mapping between states
/// and their on-node amino-acid index (which may very well represent something other
/// than the index of the enumeration element for the rotamer's "aa()").
int
FixedBBInteractionGraph::aatype_for_node_state( int, int ) const
{
	utility_exit_with_message("Unimplemented aatype_for_node_state" );
	return 0;
}

/// @details Do not call this function unless it is implemented in the derived
/// class, as should be indicated by the "aa_submatrix_energies_retrivable" method.
ObjexxFCL::FArray2D< core::PackerEnergy >
FixedBBInteractionGraph::get_aa_submatrix_energies_for_edge(
	int ,
	int ,
	int ,
	int
) const
{
	utility_exit_with_message("Unimplemented get_aa_submatrix_energies_for_edge" );
	return ObjexxFCL::FArray2D< core::PackerEnergy >();
}


} //end namespace interaction_graph
} //end namespace pack
} //end namespace core
