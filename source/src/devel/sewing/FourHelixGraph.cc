// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FourHelixGraph.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <devel/sewing/FourHelixGraph.hh>
#include <devel/sewing/Node.hh>

//Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

//Basic headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//Boost



namespace devel {
namespace sewing {

static thread_local basic::Tracer TR( "FourHelixGraph" );

FourHelixGraph::FourHelixGraph(
	utility::sql_database::sessionOP db_session_
):
ProteinAssemblyGraph(db_session_)
{}

void
FourHelixGraph::populate_graph(
	core::Real max_rmsd,
	core::Real min_clash_score
){

	std::string nodes_select =
		"SELECT\n"
		"	hgn.struct_id, hgn.bundle_id, hgn.node_id,\n"
		"	bh1.helix_id, bh1.residue_begin, bh1.residue_end,\n"
		"	bh2.helix_id, bh2.residue_begin, bh2.residue_end\n"
		"FROM helix_graph_nodes hgn\n"
		"JOIN bundle_helices bh1 ON\n"
		"	hgn.helix_id_1 = bh1.helix_id\n"
		"JOIN bundle_helices bh2 ON\n"
		"	hgn.helix_id_2 = bh2.helix_id\n"
		"ORDER BY hgn.struct_id, hgn.bundle_id;";
		
	cppdb::statement nodes_select_stmt = basic::database::safely_prepare_statement(nodes_select, db_session_);
	TR << "about to read nodes from the db" << std::endl;
	cppdb::result nodes_res = basic::database::safely_read_from_database(nodes_select_stmt);
		
	std::map<core::Size, NodeOP> nodes_map;
	//std::set<NodeOP> bundle_nodes;
	utility::vector1<NodeOP> bundle_nodes;
	core::Size prev_bundle_id=0;
	while(nodes_res.next())
	{
		protocols::features::StructureID struct_id;
		core::Size bundle_id, node_id;
		core::Size helix_1_id, helix_1_begin, helix_1_end;
		core::Size helix_2_id, helix_2_begin, helix_2_end;
		
		nodes_res >>
			struct_id >> bundle_id >> node_id >>
			helix_1_id >> helix_1_begin >> helix_1_end >>
			helix_2_id >> helix_2_begin >> helix_2_end;

		NodeOP node = new Node(
			struct_id,
			bundle_id,
			node_id,
			helix_1_id,
			helix_1_begin,
			helix_1_end,
			helix_2_id,
			helix_2_begin,
			helix_2_end);
			
		nodes_map.insert(std::make_pair(node_id, node));
		
		if(bundle_id != prev_bundle_id && prev_bundle_id != 0)
		{
			for(core::Size i=1; i<=bundle_nodes.size(); ++i)
			{
				NodeOP node_1=bundle_nodes[i];
				for(core::Size j=i+1; j<=bundle_nodes.size(); ++j)
				{
					NodeOP node_2=bundle_nodes[j];
					
					//ensure nodes share at least one helix
					if( node_1->helix_1_id() != node_2->helix_1_id() && node_1->helix_1_id() != node_2->helix_2_id() &&
						node_1->helix_2_id() != node_2->helix_1_id() && node_1->helix_2_id() != node_2->helix_2_id()) continue;
					
					//ensure nodes aren't the same helices
					if(node_1->helix_1_id() == node_2->helix_2_id() && node_1->helix_2_id() == node_2->helix_1_id()) continue;
					
					addIntraStructureEdge(node_1, node_2);
					addIntraStructureEdge(node_2, node_1);
				}
			}
			bundle_nodes.clear();
		}
		prev_bundle_id=bundle_id;
		bundle_nodes.push_back(node);
	}
	
	//record the last bundle
	for(core::Size i=1; i<=bundle_nodes.size(); ++i)
	{
		NodeOP node_1=bundle_nodes[i];
		for(core::Size j=i+1; j<=bundle_nodes.size(); ++j)
		{
			NodeOP node_2=bundle_nodes[j];
			
			//ensure nodes share at least one helix
			if( node_1->helix_1_id() != node_2->helix_1_id() && node_1->helix_1_id() != node_2->helix_2_id() &&
				node_1->helix_2_id() != node_2->helix_1_id() && node_1->helix_2_id() != node_2->helix_2_id()) continue;
			
			//ensure nodes aren't the same helices
			if(node_1->helix_1_id() == node_2->helix_2_id() && node_1->helix_2_id() == node_2->helix_1_id()) continue;
			
			addIntraStructureEdge(node_1, node_2);
			addIntraStructureEdge(node_2, node_1);
		}
	}
	
	TR << "Done making intra-bundle edges for " << nodes_map.size() << " nodes." << std::endl;
	
	std::string edges_select =
		"SELECT\n"
		"	node_id_1, node_id_2, rmsd, clash_score\n"
		"FROM node_comparisons\n"
		"WHERE\n"
		"	rmsd <= ? AND\n"
		"	clash_score >= ?;";
		
	cppdb::statement edges_select_stmt = basic::database::safely_prepare_statement(edges_select, db_session_);
	edges_select_stmt.bind(1,max_rmsd);
	edges_select_stmt.bind(2,min_clash_score);
	cppdb::result edges_res = basic::database::safely_read_from_database(edges_select_stmt);
		
	core::Size edge_counter=0;
	while(edges_res.next())
	{
		core::Size node_id_1, node_id_2;
		core::Real rmsd, clash_score;
		
		edges_res >>
			node_id_1 >> node_id_2 >> rmsd >> clash_score;
			
		NodeOP node_1 = nodes_map.find(node_id_1)->second;
		NodeOP node_2 = nodes_map.find(node_id_2)->second;
		addInterStructureEdge(node_1, node_2);
		addInterStructureEdge(node_2, node_1);
		edge_counter++;
	}
	TR << "Done making " << edge_counter << " inter-bundle edges" << std::endl;
	
}

core::Size
FourHelixGraph::get_path_size(
	core::Size num_desired_helices
){
	if(num_desired_helices % 2 != 0){
		utility_exit_with_message("An odd number of helices was given to the 4 helix graph, you can't do that!");
	}
	return ((num_desired_helices/2)*3)+1;
}

void
FourHelixGraph::tree_finder(
		EdgeList,
		utility::vector1<NodeOP>,
		core::Size /* num_helices*/,
		core::Size /*dup_counter*/,
		bool /*finish_tree*/,
		EdgeSet /*possible_inter_edges*/,
		EdgeSet /*possible_intra_edges*/
){
	utility_exit_with_message("Four helix graph is currently broken. You should fix it!");
//	if(num_helices % 2 != 0){
//		utility_exit_with_message("An odd number of helices was given to the 4 helix graph, you can't do that");
//	}
//	core::Size path_size=((num_helices/2)*3)-2;
//	
//	NodeOP last_added_node = visited_nodes[visited_nodes.size()];
//	
//	//Can only go to another node in the same bundle
//	if(visited_nodes.size() % 3 != 1){
//		std::set<NodeOP> possible_next_nodes = intra_structure_adjacency_map_[last_added_node];
//		
//		//If we are looking for the second intra bundle node then make sure it shares no helices
//		//with the originally matched pair of helices for this node
//		if(visited_nodes.size()%3==0)
//		{
//			std::set<core::Size> matched_helices;
//			matched_helices.insert(visited_nodes[visited_nodes.size()-1]->helix_1_id());
//			matched_helices.insert(visited_nodes[visited_nodes.size()-1]->helix_2_id());
//			std::set<NodeOP> temp;
//			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
//			{
//				NodeOP const & cur_node=(*it);
//				if( matched_helices.find(cur_node->helix_1_id()) == matched_helices.end() &&
//					matched_helices.find(cur_node->helix_2_id()) == matched_helices.end())
//				{
//					temp.insert(cur_node);
//				}
//			}
//			possible_next_nodes=temp;
//		}
//		
//		if(visited_nodes.size()+1 == path_size){
//			//Add the last two helices to the structure, looping through all possible next nodes
//			//and adding them here will result in two identical trees with the last nodes simply
//			//being swaps of one-another (ie pair 1-2 and pair 2-1)
//			NodeOP const & first_node = (*possible_next_nodes.begin());
//			if(!visited_nodes.contains(first_node)){
//				
//				visited_edges.push_back(Edge(last_added_node, first_node));
//				visited_nodes.push_back(first_node);
//				tree_list_.push_back(visited_edges);
//			}
//		}
//		else{
//			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
//			{
//				NodeOP const & cur_node=(*it);
//				if (visited_nodes.contains(cur_node)){
//					continue;
//				}
//				
//				visited_edges.push_back(Edge(last_added_node, cur_node));
//				visited_nodes.push_back(cur_node);
//				
//				tree_finder(visited_edges, visited_nodes, num_helices, possible_inter_edges);
//				
//				visited_nodes.pop_back();
//				visited_edges.pop_back();
//			}
//		}
//	}
//	else{
//		//Add all new edges from the last added node
//		std::set<NodeOP> const & possible_next_nodes = inter_structure_adjacency_map_[last_added_node];
//		
//		//Dont revisit a node we've been to already.
//		std::set<NodeOP> illegal_nodes;
//		illegal_nodes.insert(visited_nodes.begin(), visited_nodes.end());
//		
//		//Don't visit a pair that comes from the same bundle as any of the previously visited nodes
//		for(core::Size i=1; i<=visited_nodes.size(); ++i){
//			std::set<NodeOP> const & neighbors = intra_structure_adjacency_map_[visited_nodes[i]];
//			illegal_nodes.insert(neighbors.begin(), neighbors.end());
//		}
//		
//		for(std::set<NodeOP>::const_iterator it = possible_next_nodes.begin();
//			it != possible_next_nodes.end(); ++it){
//			possible_inter_edges.insert(Edge(last_added_node, *it));
//		}
//		
//		if(visited_nodes.size()+1 == path_size){
//			utility_exit_with_message("The path size ends on an inter-bundle match, something is wrong with helical assembly");
//		}
//		else{
//			for(EdgeSet::const_iterator it = possible_inter_edges.begin();
//				it != possible_inter_edges.end(); ++it){
//				if (illegal_nodes.find(it->second) != illegal_nodes.end()){
//					continue;
//				}
//				
//				visited_nodes.push_back(it->second);
//				visited_edges.push_back(*it);
//				
//				//Now trim nodes that were connected to the previous node, but are not attached to the one we just added.
//				//We do this to prevent putting helices on top of eachother.
//				EdgeSet possible_inter_edges_trimmed = possible_inter_edges;
//				if(visited_nodes.size() > 1){
//					
//					std::set<NodeOP> const & old_possibilities = inter_structure_adjacency_map_[it->first];//connected to previously visited node
//					std::set<NodeOP> const & new_possibilities = inter_structure_adjacency_map_[it->second];//connected to node we just added
//					
//					for(EdgeSet::const_iterator trimmed_it = possible_inter_edges.begin();
//						trimmed_it != possible_inter_edges.end(); ++trimmed_it){
//						
//						if( old_possibilities.find(trimmed_it->second) != old_possibilities.end()
//							&& new_possibilities.find(trimmed_it->second) == new_possibilities.end())
//						{
//							possible_inter_edges_trimmed.erase(*trimmed_it);
//						}
//					}
//				}
//				
//				tree_finder(visited_edges, visited_nodes, num_helices, possible_inter_edges_trimmed);
//				
//				visited_nodes.pop_back();
//				visited_edges.pop_back();
//			}
//		}
//	}
}

} //sewing namespace
} //devel namespace
