// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/Graph.hh
/// @brief  generic graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_graph_DisjointSets_hh
#define INCLUDED_utility_graph_DisjointSets_hh

// Project Headers

// Utility Headers

// C++ headers
#include <map>

#include <utility/vector1.hh>


namespace utility {
namespace graph {

struct DS_Node {
	platform::Size index;
	platform::Size parent;
	platform::Size rank;
};

class DisjointSets
{

public:

	DisjointSets();

	///  @brief ctor for class if the number of nodes is known upfront. Fastest.
	DisjointSets( platform::Size n_nodes );

	platform::Size n_nodes() const;

	/// @brief add a new node -- as implemented, O(N) operation
	/// use the DS(platform::Size) constructor for better speed.
	void ds_make_set();

	/// @brief return the representative for a node
	platform::Size ds_find( platform::Size node_id ) const;

	/// @brief make it so that two nodes end up in the same set
	void ds_union( platform::Size node1, platform::Size node2 );

	/// @brief DS_Node read access
	DS_Node const &
	node( platform::Size node_id ) const {
		return nodes_[ node_id ];
	}

	/// @brief count the number of disjoint sets. O(N)
	platform::Size n_disjoint_sets() const;

	/// @brief count the size of each disjoint set. O(N)
	utility::vector1< platform::Size >
	disjoint_set_sizes() const;

	/// @brief return the nodes in the set containing the specified node.
	utility::vector1< platform::Size >
	nodes_in_set( platform::Size node_id ) const;

	/// @brief return a map from the representative node of each set to
	///  the list of nodes in their sets
	std::map< platform::Size, utility::vector1< platform::Size > >
	sets() const;

private:

	mutable utility::vector1< DS_Node > nodes_;

};


}
}

#endif
