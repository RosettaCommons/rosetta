// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FourHelixRepeatGraph.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_FourHelixRepeatGraph_HH
#define INCLUDED_devel_sewing_FourHelixRepeatGraph_HH

#include <devel/sewing/FourHelixGraph.hh>

#include <core/types.hh>

#include <utility/vector1.fwd.hh>

namespace devel {
namespace sewing {

class FourHelixRepeatGraph : public FourHelixGraph {
public:
	
	FourHelixRepeatGraph(
		utility::sql_database::sessionOP db_session
	);
	
	core::Size
	get_path_size(
		core::Size num_desired_helices
	);
	
	void
	find_paths(
		EdgeList visited_edges,
		utility::vector1<NodeOP> visited_nodes,
		core::Size path_size
	);
	
	void
	tree_finder(
		EdgeList visited_edges,
		utility::vector1<NodeOP> visited_nodes,
		core::Size path_size,
		core::Size dup_counter,
		bool finish_tree,
		EdgeSet possible_inter_edges,
		EdgeSet possible_intra_edges
	);
	
private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size num_cycles_to_print_;
};

} //sewing namespace
} //devel namespace

#endif
