// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ThreeHelixGraph.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_ThreeHelixGraph_hh
#define INCLUDED_ThreeHelixGraph_hh

//Unit headers
#include <devel/sewing/ThreeHelixGraph.fwd.hh>
#include <devel/sewing/ProteinAssemblyGraph.hh>
#include <devel/sewing/Node.hh>

//Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//Boost


//C++ headers
#include <map>
#include <set>
#include <string>

namespace devel {
namespace sewing {

class ThreeHelixGraph : public ProteinAssemblyGraph
{
public:

	ThreeHelixGraph(
		utility::sql_database::sessionOP db_session
	);

	virtual core::Size
	get_path_size(
		core::Size num_desired_helices
	);

	virtual void
	tree_finder(
		EdgeList visited_edges,
		utility::vector1<NodeOP> visited_nodes,
		core::Size num_helices,
		core::Size dup_counter,
		bool finish_tree,
		EdgeSet possible_inter_edges,
		EdgeSet possible_intra_edges
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
