// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixAssemblyGraph.hh
///
/// @brief
/// @author Tim Jacobs



#ifndef INCLUDED_HelixAssemblyGraph_hh
#define INCLUDED_HelixAssemblyGraph_hh

//Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//Utility headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++ headers
#include <map>
#include <set>
#include <string>

class Node
{
	
public:
		
	Node(core::Size pair_id);
	
	core::Size pair_id() const;
	std::string print() const;
	
	bool operator<(const Node & rhs) const;
	bool operator==(const Node & rhs) const;
	
private:
	core::Size pair_id_;
	
};

class HelixAssemblyGraph
{
public:

	typedef std::map< Node, utility::vector1<Node> > AdjacencyMap;
	typedef utility::vector1< utility::vector1<Node> > PathList;

	HelixAssemblyGraph(utility::sql_database::sessionOP db_session, bool simple_print);

	void addEdge(Node node1, Node node2);
	
	void addIntraBundleEdge(Node node1, Node node2);
	void addInterBundleEdge(Node node1, Node node2);
	
	void printPathSimple(utility::vector1<Node> visited);
	
	void printPath(utility::vector1<Node> visited);
	
	void combinePoses(core::pose::Pose & pose1, core::pose::Pose const & pose2);
	
	void depthFirst(utility::vector1<Node> visited, size_t path_size);

	core::Real bb_score(core::pose::Pose & pose, core::Size new_helix_begin, core::Size new_helix_end);

	std::set<Node> nodes();
	
	core::Size total_inter_bundle_edges();
	
private:
	AdjacencyMap inter_bundle_adjacency_map_;
	AdjacencyMap intra_bundle_adjacency_map_;
	std::set<Node> nodes_;
	utility::sql_database::sessionOP db_session_;
	PathList path_list_;
	bool simple_print_;
	core::Size inter_bundle_edge_counter_;
};

#endif
