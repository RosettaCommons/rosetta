// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FourHelixRepeatGraph.cc
///
/// @brief
/// @author Tim Jacobs

//Unit Headers
#include <devel/sewing/FourHelixRepeatGraph.hh>

//Core Headers
#include <core/types.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

//Basic Headers
#include <basic/Tracer.hh>

namespace devel {
namespace sewing {

static thread_local basic::Tracer TR( "FourHelixRepeatGraph" );

FourHelixRepeatGraph::FourHelixRepeatGraph(
	utility::sql_database::sessionOP db_session
):
FourHelixGraph(db_session)
{}

core::Size
FourHelixRepeatGraph::get_path_size(
	core::Size num_desired_helices
){
	if(num_desired_helices % 2 != 0){
		utility_exit_with_message("An odd number of helices was given to the 4 helix repeat graph, you can't do that!");
	}
	return ((num_desired_helices/2)*3)+1;
}

void
FourHelixRepeatGraph::find_paths(
	EdgeList visited_edges,
	utility::vector1<NodeOP> visited_nodes,
	core::Size path_size
){
	NodeOP last_added_node = visited_nodes[visited_nodes.size()];
	
	//Can only go to another node in the same bundle
	if(visited_nodes.size() % 3 != 1){
		std::set<NodeOP> possible_next_nodes = intra_structure_adjacency_map_[last_added_node];
		
		//If we are looking for the second intra bundle node then make sure it shares no helices
		//with the originally matched pair of helices for this node
		if(visited_nodes.size()%3==0)
		{
			std::set<core::Size> matched_helices;
			matched_helices.insert(visited_nodes[visited_nodes.size()-1]->helix_1_id());
			matched_helices.insert(visited_nodes[visited_nodes.size()-1]->helix_2_id());
			std::set<NodeOP> temp;
			NodeOP cur_node;
			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
			{
				NodeOP cur_node=(*it);
				if( matched_helices.find(cur_node->helix_1_id()) == matched_helices.end() &&
					matched_helices.find(cur_node->helix_2_id()) == matched_helices.end())
				{
					temp.insert(cur_node);
				}
			}
			possible_next_nodes=temp;
		}
		
		if(visited_nodes.size()+1 == path_size)
		{
			NodeOP cur_node;
			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
			{
				cur_node=(*it);
				if(visited_nodes[1]==cur_node){
					//don't add the last edge, it's only there to make sure we have a full repeat
					visited_edges.push_back(Edge(last_added_node, cur_node));
					visited_nodes.push_back(cur_node);
					tree_list_.push_back(visited_edges);
				}
			}
		}
		else{
			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
			{
				NodeOP cur_node=(*it);
				if (visited_nodes.contains(cur_node)){
					continue;
				}
				
				visited_edges.push_back(Edge(last_added_node, cur_node));
				visited_nodes.push_back(cur_node);
				
				find_paths(visited_edges, visited_nodes, path_size);
				
				visited_nodes.pop_back();
				visited_edges.pop_back();
			}
		}
	}
	else{
		//Add all new edges from the last added node
		std::set<NodeOP> const & possible_next_nodes = inter_structure_adjacency_map_[last_added_node];
		
		//Dont revisit a node we've been to already.
		std::set<NodeOP> illegal_nodes;
		illegal_nodes.insert(visited_nodes.begin(), visited_nodes.end());
		
		//Don't visit a pair that comes from the same bundle as any of the previously visited nodes,
		//unless this is the step to cycle back and cause the repeat
		if(visited_nodes.size()+3 != path_size){
			for(core::Size i=1; i<=visited_nodes.size(); ++i){
				std::set<NodeOP> const & neighbors = intra_structure_adjacency_map_[visited_nodes[i]];
				illegal_nodes.insert(neighbors.begin(), neighbors.end());
			}
		}
		
		if(visited_nodes.size()+1 == path_size){
			utility_exit_with_message("The path size ends on an inter-bundle match, something is wrong with helical assembly");
		}
		else{
			NodeOP cur_node;
			for(std::set<NodeOP>::const_iterator it=possible_next_nodes.begin(); it !=possible_next_nodes.end(); ++it)
			{
				cur_node=(*it);
				if (illegal_nodes.find(cur_node) != illegal_nodes.end())
				{
					continue;
				}
				
				visited_edges.push_back(Edge(last_added_node, cur_node));
				visited_nodes.push_back(cur_node);
				
				find_paths(visited_edges, visited_nodes, path_size);
				
				visited_nodes.pop_back();
				visited_edges.pop_back();
			}
		}
	}
}
	
void
FourHelixRepeatGraph::tree_finder(
	EdgeList visited_edges,
	utility::vector1<NodeOP> visited_nodes,
	core::Size path_size,
	core::Size /*dup_counter*/,
	bool /*finish_tree*/,
	EdgeSet /*possible_inter_edges*/,
	EdgeSet /*possible_intra_edges*/
){
	find_paths(visited_edges, visited_nodes, path_size);
}

} //sewing namespace
} //devel namespace
