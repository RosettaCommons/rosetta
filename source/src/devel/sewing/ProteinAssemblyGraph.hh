// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProteinAssemblyGraph.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_ProteinAssemblyGraph_hh
#define INCLUDED_devel_sewing_ProteinAssemblyGraph_hh

//Unit headers
#include <devel/sewing/ProteinAssemblyGraph.fwd.hh>
#include <devel/sewing/Node.hh>

//Protocols
#include <protocols/features/FeaturesReporter.fwd.hh>

//Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

//Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

//C++ headers
#include <map>
#include <set>
#include <string>

namespace devel {
namespace sewing {

class ProteinAssemblyGraph : public utility::pointer::ReferenceCount
{
public:

	typedef std::map< NodeOP, std::set<NodeOP> > AdjacencyMap;
	typedef utility::vector1< std::set<NodeOP> > PathList;
	typedef std::pair<NodeOP, NodeOP> Edge;
	typedef utility::vector1< Edge > EdgeList;
	typedef std::set< Edge > EdgeSet;

	ProteinAssemblyGraph(utility::sql_database::sessionOP db_session);

	void addIntraStructureEdge(NodeOP node1, NodeOP node2);
	void addInterStructureEdge(NodeOP node1, NodeOP node2);

	void
	print_tree_simple(
		EdgeList visited
	);

	//Given a list of edges, create structures
	void
	print_tree(
		EdgeList visited,
		core::Size tree_id
	);

	utility::vector1<EdgeList>
	find_all_trees(
		core::Size num_helices
	);

	virtual core::Size
	get_path_size(
		core::Size num_desired_helices
	) = 0;

	virtual void
	tree_finder(
		EdgeList visited_edges,
		utility::vector1<NodeOP> visited_nodes,
		core::Size path_size,
		core::Size dup_counter,
		bool finish_tree,
		EdgeSet possible_inter_edges,
		EdgeSet possible_intra_edges
	) = 0;

	virtual void
	populate_graph(
		core::Real max_rmsd,
		core::Real min_clash_score
	)=0;

	core::Real
	bb_score(
		core::pose::Pose & pose,
		core::Size new_residues_begin,//change to take a set of residue numbers
		core::Size new_residues_end
	);


	virtual core::pose::Pose
	print_first(
		EdgeList const & visited_edges,
		std::map<NodeOP, core::pose::Pose> const & poses
	);

	void
	dump_resfile(
		std::map<core::Size, utility::vector1<char> > sequence_map,
		std::string filename
	);

	std::set<NodeOP> nodes();

	core::Size total_inter_structure_edges();

	utility::vector1<EdgeList> tree_list();

protected:
	AdjacencyMap inter_structure_adjacency_map_;
	AdjacencyMap intra_structure_adjacency_map_;
	std::set<NodeOP> nodes_;
	utility::sql_database::sessionOP db_session_;
	PathList path_list_;
	utility::vector1<EdgeList> tree_list_;

	//cache poses to increase speed of assembly
	std::map<NodeOP, core::pose::Pose> poses_;

	core::Size inter_structure_edge_counter_;

	core::scoring::ScoreFunctionOP scorefxn_;
};

} //sewing namespace
} //devel namespace

#endif
