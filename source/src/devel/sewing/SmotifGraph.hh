// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SmotifGraph.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_SmotifGraph_HH
#define INCLUDED_devel_sewing_SmotifGraph_HH

//Unit Headers
#include <devel/sewing/SmotifGraph.fwd.hh>
#include <devel/sewing/ProteinAssemblyGraph.hh>
#include <devel/sewing/Node.hh>

//Utility Headers
#include <utility/vector1.fwd.hh>

namespace devel {
namespace sewing {

class SmotifGraph : public ProteinAssemblyGraph
{

public:

	SmotifGraph(
		utility::sql_database::sessionOP db_session
	);
	
	void
	print_tree(
		EdgeList visited,
		core::Size tree_id
	);

	virtual void
	tree_finder(
		EdgeList visited_edges,
		utility::vector1<Node> visited_nodes,
		core::Size num_helices,
		EdgeSet possible_list
	);
	
	void
	populate_graph(
		core::Real max_rmsd,
		core::Real min_clash_score
	);
	
};

} //sewing namespace
} //devel namespace

#endif
