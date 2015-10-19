// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SmotifGraph.cc
///
/// @brief
/// @author Tim Jacobs

#include <devel/sewing/SmotifGraph.hh>

//Basic headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

SmotifGraph::SmotifGraph(
	utility::sql_database::sessionOP db_session_
):
ProteinAssemblyGraph(db_session_)
{}

void
SmotifGraph::populate_graph(
	core::Real max_rmsd,
	core::Real min_clash_score
){
	std::string inter_structure_select_string =
		"SELECT\n"
		"	hgn1.struct_id, hgn1.bundle_id, hgn1.node_id,\n"
		"	bh1.helix_id, bh1.residue_begin, bh1.residue_end,\n"
		"	bh2.helix_id, bh2.residue_begin, bh2.residue_end,\n"
		"	hgn2.struct_id, hgn2.bundle_id, hgn2.node_id,\n"
		"	bh3.helix_id, bh3.residue_begin, bh3.residue_end,\n"
		"	bh4.helix_id, bh4.residue_begin, bh4.residue_end,\n"
		"	c.rmsd\n"
		"FROM node_comparisons c\n"
		"JOIN helix_graph_nodes hgn1 ON\n"
		"	c.node_id_1 = hgn1.node_id\n"
		"JOIN bundle_helices bh1 ON\n"
		"	hgn1.helix_id_1 = bh1.helix_id\n"
		"JOIN bundle_helices bh2 ON\n"
		"	hgn1.helix_id_2 = bh2.helix_id\n"
		"JOIN helix_graph_nodes hgn2 ON\n"
		"	c.node_id_2 = hgn2.node_id\n"
		"JOIN bundle_helices bh3 ON\n"
		"	hgn2.helix_id_1 = bh3.helix_id\n"
		"JOIN bundle_helices bh4 ON\n"
		"	hgn2.helix_id_2 = bh4.helix_id\n"
		"WHERE\n"
		"	rmsd <= "+utility::to_string(max_rmsd)+"AND\n"
		"	clash_score > "+utility::to_string(min_clash_score)+";";
	
	cppdb::statement inter_structure_select_stmt = basic::database::safely_prepare_statement(inter_structure_select_string, db_session_);
	cppdb::result res1 = basic::database::safely_read_from_database(inter_structure_select_stmt);
	
	while(res1.next()){
		protocols::features::StructureID struct_id_1;
		core::Size bundle_id_1, node_id_1;
		core::Size pair_1_helix_1_id, pair_1_helix_1_begin, pair_1_helix_1_end;
		core::Size pair_1_helix_2_id, pair_1_helix_2_begin, pair_1_helix_2_end;
		
		protocols::features::StructureID struct_id_2;
		core::Size bundle_id_2, node_id_2;
		core::Size pair_2_helix_1_id, pair_2_helix_1_begin, pair_2_helix_1_end;
		core::Size pair_2_helix_2_id, pair_2_helix_2_begin, pair_2_helix_2_end;
		core::Real rmsd;
		
		res1 >>
			struct_id_1 >> bundle_id_1 >> node_id_1 >>
			pair_1_helix_1_id >> pair_1_helix_1_begin >> pair_1_helix_1_end >>
			pair_1_helix_2_id >> pair_1_helix_2_begin >> pair_1_helix_2_end >>
			struct_id_2 >> bundle_id_2 >> node_id_2 >>
			pair_2_helix_1_id >> pair_2_helix_1_begin >> pair_2_helix_1_end >>
			pair_2_helix_2_id >> pair_2_helix_2_begin >> pair_2_helix_2_end >>
			rmsd;
		
		Node node1(
			struct_id_1,
			bundle_id_1,
			node_id_1,
			pair_1_helix_1_id,
			pair_1_helix_1_begin,
			pair_1_helix_1_end,
			pair_1_helix_2_id,
			pair_1_helix_2_begin,
			pair_1_helix_2_end);
		
		Node node2(
			struct_id_2,
			bundle_id_2,
			node_id_2,
			pair_2_helix_1_id,
			pair_2_helix_1_begin,
			pair_2_helix_1_end,
			pair_2_helix_2_id,
			pair_2_helix_2_begin,
			pair_2_helix_2_end);
		
		addInterStructureEdge(node1, node2);
		addInterStructureEdge(node2, node1);
	}
	
	std::string intra_bundle_select_string =
		"SELECT\n"
		"	hgn1.struct_id, hgn1.bundle_id, hgn1.node_id,\n"
		"	bh1.helix_id, bh1.residue_begin, bh1.residue_end,\n"
		"	bh2.helix_id, bh2.residue_begin, bh2.residue_end,\n"
		"	hgn2.struct_id, hgn2.bundle_id, hgn2.node_id,\n"
		"	bh3.helix_id, bh3.residue_begin, bh3.residue_end,\n"
		"	bh4.helix_id, bh4.residue_begin, bh4.residue_end\n"
		"FROM helix_graph_nodes hgn1\n"
		"JOIN bundle_helices bh1 ON\n"
		"	hgn1.helix_id_1 = bh1.helix_id\n"
		"JOIN bundle_helices bh2 ON\n"
		"	hgn1.helix_id_2 = bh2.helix_id\n"
		"JOIN helix_graph_nodes hgn2 ON\n"
		"	hgn1.bundle_id=hgn2.bundle_id AND\n"
		"	hgn1.node_id < hgn2.node_id AND\n"
		"	(hgn1.helix_id_1 <> hgn2.helix_id_2 OR hgn1.helix_id_2 <> hgn2.helix_id_1)"
		"JOIN bundle_helices bh3 ON\n"
		"	hgn2.helix_id_1 = bh3.helix_id\n"
		"JOIN bundle_helices bh4 ON\n"
		"	hgn2.helix_id_2 = bh4.helix_id\n";
	
	cppdb::statement intra_bundle_select_stmt = basic::database::safely_prepare_statement(intra_bundle_select_string, db_session_);
	cppdb::result res2 = basic::database::safely_read_from_database(intra_bundle_select_stmt);
	
	while(res2.next()){
		protocols::features::StructureID struct_id_1;
		core::Size bundle_id_1, node_id_1;
		core::Size pair_1_helix_1_id, pair_1_helix_1_begin, pair_1_helix_1_end;
		core::Size pair_1_helix_2_id, pair_1_helix_2_begin, pair_1_helix_2_end;
		
		protocols::features::StructureID struct_id_2;
		core::Size bundle_id_2, node_id_2;
		core::Size pair_2_helix_1_id, pair_2_helix_1_begin, pair_2_helix_1_end;
		core::Size pair_2_helix_2_id, pair_2_helix_2_begin, pair_2_helix_2_end;
		
		res2 >>
		struct_id_1 >> bundle_id_1 >> node_id_1 >>
		pair_1_helix_1_id >> pair_1_helix_1_begin >> pair_1_helix_1_end >>
		pair_1_helix_2_id >> pair_1_helix_2_begin >> pair_1_helix_2_end >>
		struct_id_2 >> bundle_id_2 >> node_id_2 >>
		pair_2_helix_1_id >> pair_2_helix_1_begin >> pair_2_helix_1_end >>
		pair_2_helix_2_id >> pair_2_helix_2_begin >> pair_2_helix_2_end;
		
		Node node1(
			struct_id_1,
			bundle_id_1,
			node_id_1,
			pair_1_helix_1_id,
			pair_1_helix_1_begin,
			pair_1_helix_1_end,
			pair_1_helix_2_id,
			pair_1_helix_2_begin,
			pair_1_helix_2_end);
		
		Node node2(
			struct_id_2,
			bundle_id_2,
			node_id_2,
			pair_2_helix_1_id,
			pair_2_helix_1_begin,
			pair_2_helix_1_end,
			pair_2_helix_2_id,
			pair_2_helix_2_begin,
			pair_2_helix_2_end);
		
		addIntraStructureEdge(node1, node2);
		addIntraStructureEdge(node2, node1);
	}
}

void
SmotifGraph::tree_finder(
	EdgeList visited_edges,
	utility::vector1<Node> visited_nodes,
	core::Size num_helices,
	EdgeSet possible_edges
){
	core::Size path_size=(num_helices*2)-3;

	Node last_added_node = visited_nodes[visited_nodes.size()];

	//Can only go to another node in the same bundle
	if(visited_nodes.size() % 2 == 0){
		utility::vector1<Node> possible_next_nodes = intra_structure_adjacency_map_[last_added_node];

		if(visited_nodes.size()+1 == path_size){
			if(!visited_nodes.contains(possible_next_nodes[1])){

				visited_edges.push_back(Edge(last_added_node, possible_next_nodes[1]));
				visited_nodes.push_back(possible_next_nodes[1]);
				tree_list_.push_back(visited_edges);
			}
		}
		else{
			for(core::Size i=1; i<=possible_next_nodes.size(); ++i){
				if (visited_nodes.contains(possible_next_nodes[i])){
					continue;
				}

				visited_edges.push_back(Edge(last_added_node, possible_next_nodes[i]));
				visited_nodes.push_back(possible_next_nodes[i]);

				tree_finder(visited_edges, visited_nodes, num_helices, possible_edges);

				visited_nodes.pop_back();
				visited_edges.pop_back();
			}
		}
	}
	else{
		//Add all new edges from the last added node
		std::set<Node> possible_next_nodes;
		possible_next_nodes.insert(inter_structure_adjacency_map_[last_added_node].begin(), inter_structure_adjacency_map_[last_added_node].end());

		//Dont revisit a node we've been to already. 
		std::set<Node> illegal_nodes;
		illegal_nodes.insert(visited_nodes.begin(), visited_nodes.end());
		
		//Don't visit a pair that comes from the same bundle as any of the previously visited nodes
		for(core::Size i=1; i<=visited_nodes.size(); ++i){
			utility::vector1<Node> neighbors = intra_structure_adjacency_map_[visited_nodes[i]];
			illegal_nodes.insert(neighbors.begin(), neighbors.end());
		}

		for(std::set<Node>::const_iterator it = possible_next_nodes.begin();
			it != possible_next_nodes.end(); ++it){
			possible_edges.insert(Edge(last_added_node, *it));
		}

		if(visited_nodes.size()+1 == path_size){
			utility_exit_with_message("The path size ends on an inter-bundle match, something is wrong with helical assembly");
		}
		else{
			for(EdgeSet::const_iterator it = possible_edges.begin();
				it != possible_edges.end(); ++it){
				if (illegal_nodes.find(it->second) != illegal_nodes.end()){
					continue;
				}

				visited_nodes.push_back(it->second);
				visited_edges.push_back(*it);

				//Now trim nodes that were connected to the previous node, but are not attached to the one we just added.
				//We do this to prevent putting helices on top of eachother.
				EdgeSet possible_edges_trimmed = possible_edges;
				if(visited_nodes.size() > 1){

					utility::vector1<Node> old_possibilities = inter_structure_adjacency_map_[it->first];//connected to previously visited node
					utility::vector1<Node> new_possibilities = inter_structure_adjacency_map_[it->second];//connected to node we just added
					
					for(EdgeSet::const_iterator trimmed_it = possible_edges.begin();
						trimmed_it != possible_edges.end(); ++trimmed_it){
						
						if(old_possibilities.contains(trimmed_it->second) && !new_possibilities.contains(trimmed_it->second)){
							possible_edges_trimmed.erase(*trimmed_it);
						}
					}
				}

				tree_finder(visited_edges, visited_nodes, num_helices, possible_edges_trimmed);

				visited_nodes.pop_back();
				visited_edges.pop_back();
			}
		}
	}
}

void
SmotifGraph::print_tree(
	EdgeList visited,
	core::Size tree_id
){

	std::string output_name = "tree_" + utility::to_string(tree_id) + "_";
	core::pose::Pose master_pose;

	std::string select_pair_residues =
	"SELECT r.seqpos AS resnum\n"
	"FROM protein_residue_conformation r\n"
	"WHERE\n"
	"	r.struct_id = ? AND"
	"	(r.seqpos BETWEEN ? AND ? OR\n"
	"	r.seqpos BETWEEN ? AND ?)\n";

	cppdb::statement select_pair_residues_stmt =
		basic::database::safely_prepare_statement(select_pair_residues, db_session_);

	protocols::features::ProteinSilentReportOP protein_silent_report = new protocols::features::ProteinSilentReport();

	std::map<std::string, core::pose::Pose> output_poses;
	std::map<Node, core::pose::Pose> poses;

	//retrieve all the structural info we need for each bundle pair in the given tree
	for(core::Size i=1; i<=visited.size(); ++i){

		Node cur_node;
		for(core::Size j=1; j<=2; ++j){

			if(j==1){
				cur_node=visited[i].first;
			}
			else{
				cur_node=visited[i].second;
			}

			if(poses.find(cur_node) == poses.end()){//new node

				std::set<core::Size> resnums;
				core::pose::Pose pose;

				select_pair_residues_stmt.bind(1, cur_node.struct_id());
				select_pair_residues_stmt.bind(2, cur_node.helix_1_begin());
				select_pair_residues_stmt.bind(3, cur_node.helix_1_end());
				select_pair_residues_stmt.bind(4, cur_node.helix_2_begin());
				select_pair_residues_stmt.bind(5, cur_node.helix_2_end());
				cppdb::result res = basic::database::safely_read_from_database(select_pair_residues_stmt);
				while(res.next()){
					core::Size resnum;
					res >> resnum;
					resnums.insert(resnum);
				}

				protein_silent_report->load_pose(db_session_, cur_node.struct_id(), resnums, pose);
				poses[cur_node]=pose;
				
			}
		}
	}

	//Clear tree memory from previous assemblies
	for(core::Size i=1; i<=visited.size(); ++i){
		visited[i].first.helix_1_positions_.clear();
		visited[i].second.helix_1_positions_.clear();
		visited[i].first.helix_2_positions_.clear();
		visited[i].second.helix_2_positions_.clear();
	}

	//Assemble into the final bundle
	std::map<core::Size, utility::vector1<char> > sequence_map;//for resfile
	for(core::Size i=1; i<=visited.size(); ++i){

		core::Size helix_size = visited[i].first.helix_1_end()-visited[i].first.helix_1_begin()+1;

		if(i==1){
			std::string frag_name = utility::to_string(tree_id) + "_" + utility::to_string(i) + "_" + visited[i].first.print() + ".pdb";
			output_poses.insert(std::make_pair(frag_name, poses[visited[i].first]));
			
			//record residues for first helix pair
			core::Size position_counter=0;
			Node cur_node = visited[i].first;
			if(cur_node.helix_2_begin() > cur_node.helix_1_begin()){
				for(core::Size j=1; j<=helix_size; ++j){
					position_counter++;
					visited[i].first.helix_1_positions_[j]=position_counter;
				}
				for(core::Size j=1; j<=helix_size; ++j){
					position_counter++;
					visited[i].first.helix_2_positions_[j]=position_counter;
				}
			}
			else{//helices will be assembled in the opposite sequence order they appear in the native
				for(core::Size j=1; j<=helix_size; ++j){
					position_counter++;
					visited[i].first.helix_2_positions_[j]=position_counter;
				}
				for(core::Size j=1; j<=helix_size; ++j){
					position_counter++;
					visited[i].first.helix_1_positions_[j]=position_counter;
				}
			}
			
			for(std::map<core::Size,core::Size>::const_iterator it=visited[i].first.helix_1_positions_.begin();
				it!=visited[i].first.helix_1_positions_.end(); ++it){
				if(visited[i].first.reversed()){
					sequence_map[it->second].push_back(poses[visited[i].first].residue(it->first+helix_size).name1());
				}
				else{
					sequence_map[it->second].push_back(poses[visited[i].first].residue(it->first).name1());
				}
			}
			for(std::map<core::Size,core::Size>::const_iterator it=visited[i].first.helix_2_positions_.begin();
				it!=visited[i].first.helix_2_positions_.end(); ++it){
				if(visited[i].first.reversed()){
					sequence_map[it->second].push_back(poses[visited[i].first].residue(it->first).name1());
				}
				else{
					sequence_map[it->second].push_back(poses[visited[i].first].residue(it->first+helix_size).name1());
				}
			}
		}
		else{//ug, Node lists should really be NodeOP lists
			visited[i].first = visited[i-1].second;
		}

		//Pairs come from the same bundle, only superimpose the matching helix
		if(visited[i].first.bundle_id() == visited[i].second.bundle_id()){

			core::Size first_pair_offset;
			core::Size second_pair_offset;
			core::Size new_helix_start;
			//all helices should be the same size

			//Find which helix matches between the two pairs
			if(visited[i].first.helix_1_id()== visited[i].second.helix_1_id()){
				first_pair_offset=(visited[i].first.helix_1_id()< visited[i].first.helix_2_id()) ? 0 : helix_size;
				second_pair_offset=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? 0 : helix_size;
				new_helix_start=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? helix_size+1 : 1;
				
				visited[i].second.helix_1_positions_= visited[i].first.helix_1_positions_;
				core::Size position_start = master_pose.total_residue();
				for(core::Size j=1; j<=helix_size; ++j){
					visited[i].second.helix_2_positions_[j]=position_start+j;
				}
			}
			else if(visited[i].first.helix_1_id()== visited[i].second.helix_2_id()){
				first_pair_offset=(visited[i].first.helix_1_id()< visited[i].first.helix_2_id()) ? 0 : helix_size;
				second_pair_offset=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? helix_size : 0;
				new_helix_start=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? 1 : helix_size+1;
				
				visited[i].second.helix_2_positions_= visited[i].first.helix_1_positions_;
				core::Size position_start = master_pose.total_residue();
				for(core::Size j=1; j<=helix_size; ++j){
					visited[i].second.helix_1_positions_[j]=position_start+j;
				}
			}
			else if(visited[i].first.helix_2_id()== visited[i].second.helix_1_id()){
				first_pair_offset=(visited[i].first.helix_1_id()< visited[i].first.helix_2_id()) ? helix_size : 0;
				second_pair_offset=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? 0 : helix_size;
				new_helix_start=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? helix_size+1 : 1;
				
				visited[i].second.helix_1_positions_= visited[i].first.helix_2_positions_;
				core::Size position_start = master_pose.total_residue();
				for(core::Size j=1; j<=helix_size; ++j){
					visited[i].second.helix_2_positions_[j]=position_start+j;
				}
			}
			else if(visited[i].first.helix_2_id()== visited[i].second.helix_2_id()){
				first_pair_offset=(visited[i].first.helix_1_id()< visited[i].first.helix_2_id()) ? helix_size : 0;
				second_pair_offset=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? helix_size : 0;;
				new_helix_start=(visited[i].second.helix_1_id()< visited[i].second.helix_2_id()) ? 1 : helix_size+1;
				
				visited[i].second.helix_2_positions_= visited[i].first.helix_2_positions_;
				core::Size position_start = master_pose.total_residue();
				for(core::Size j=1; j<=helix_size; ++j){
					visited[i].second.helix_1_positions_[j]=position_start+j;
				}
			}
			else{
				utility_exit_with_message("No helix ids matched. something has gone wrong with helical assembly!");
				//Gets rid of some warnings
				new_helix_start=0;
				first_pair_offset=0;
				second_pair_offset=0;
			}

			core::id::AtomID_Map< core::id::AtomID > atom_map;
			atom_map.clear();
			core::pose::initialize_atomid_map( atom_map, poses[visited[i].second], core::id::BOGUS_ATOM_ID );

			for(core::Size j=1; j<= helix_size; j++){
				core::id::AtomID const id1( poses[visited[i].first].residue(j+first_pair_offset).atom_index("CA"), j+first_pair_offset);
				core::id::AtomID const id2( poses[visited[i].second].residue(j+second_pair_offset).atom_index("CA"), j+second_pair_offset);
				atom_map[ id2 ] = id1;

				core::id::AtomID const id3( poses[visited[i].first].residue(j+first_pair_offset).atom_index("C"), j+first_pair_offset);
				core::id::AtomID const id4( poses[visited[i].second].residue(j+second_pair_offset).atom_index("C"), j+second_pair_offset);
				atom_map[ id4 ] = id3;

				core::id::AtomID const id5( poses[visited[i].first].residue(j+first_pair_offset).atom_index("N"), j+first_pair_offset);
				core::id::AtomID const id6( poses[visited[i].second].residue(j+second_pair_offset).atom_index("N"), j+second_pair_offset);
				atom_map[ id6 ] = id5;

				core::id::AtomID const id7( poses[visited[i].first].residue(j+first_pair_offset).atom_index("O"), j+first_pair_offset);
				core::id::AtomID const id8( poses[visited[i].second].residue(j+second_pair_offset).atom_index("O"), j+second_pair_offset);
				atom_map[ id8 ] = id7;
			}
			core::scoring::superimpose_pose(poses[visited[i].second], poses[visited[i].first]/*const*/, atom_map);

			core::Size pose_size_before_addition = master_pose.total_residue();
			core::pose::append_subpose_to_pose(master_pose, poses[visited[i].second], new_helix_start, new_helix_start+helix_size-1, false);
			output_name+=utility::to_string(visited[i].first.pair_id()) + "_";
			
//			core::Real clash_score = bb_score(master_pose, pose_size_before_addition+1, master_pose.total_residue());
//			//		TR << "Backbone-backbone score for addition of pair " << tag_1 << ": " << clash_score << std::endl;
//			if(clash_score > 5){
//				return;
//			}
			
		}

		//Pairs come from different bundles, superimpose both helices
		else{
			visited[i].second.helix_1_positions_= visited[i].first.helix_1_positions_;
			visited[i].second.helix_2_positions_= visited[i].first.helix_2_positions_;
		
			core::id::AtomID_Map< core::id::AtomID > atom_map;
			atom_map.clear();
			core::pose::initialize_atomid_map( atom_map, poses[visited[i].second], core::id::BOGUS_ATOM_ID );

			utility::vector1<core::id::AtomID> helix_1_atoms;
			utility::vector1<core::id::AtomID> helix_2_atoms;
			
			core::Size pair_1_helix_1_offset;
			core::Size pair_1_helix_2_offset;
			core::Size pair_2_helix_1_offset;
			core::Size pair_2_helix_2_offset;
			if(visited[i].first.helix_1_begin() < visited[i].first.helix_2_begin()){
				pair_1_helix_1_offset = 0;
				pair_1_helix_2_offset = helix_size;
			}
			else{
				pair_1_helix_1_offset = helix_size;
				pair_1_helix_2_offset = 0;
			}
			
			if(visited[i].second.helix_1_begin() < visited[i].second.helix_2_begin()){
				pair_2_helix_1_offset = 0;
				pair_2_helix_2_offset = helix_size;
			}
			else{
				pair_2_helix_1_offset = helix_size;
				pair_2_helix_2_offset = 0;
			}
			
			for(core::Size j=1; j<=helix_size; j++){
				core::id::AtomID const id1( poses[visited[i].first].residue(j+pair_1_helix_1_offset).atom_index("CA"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id1);
				
				core::id::AtomID const id2( poses[visited[i].second].residue(j+pair_2_helix_1_offset).atom_index("CA"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id2);
				atom_map[ id2 ] = id1;
				
				core::id::AtomID const id3( poses[visited[i].first].residue(j+pair_1_helix_1_offset).atom_index("C"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id3);
				core::id::AtomID const id4( poses[visited[i].second].residue(j+pair_2_helix_1_offset).atom_index("C"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id4);
				atom_map[ id4 ] = id3;
				
				core::id::AtomID const id5( poses[visited[i].first].residue(j+pair_1_helix_1_offset).atom_index("N"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id5);
				core::id::AtomID const id6( poses[visited[i].second].residue(j+pair_2_helix_1_offset).atom_index("N"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id6);
				atom_map[ id6 ] = id5;
				
				core::id::AtomID const id7( poses[visited[i].first].residue(j+pair_1_helix_1_offset).atom_index("O"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id7);
				core::id::AtomID const id8( poses[visited[i].second].residue(j+pair_2_helix_1_offset).atom_index("O"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id8);
				atom_map[ id8 ] = id7;
			}
			
			for(core::Size j=1; j<=helix_size; j++){
				core::id::AtomID const id1( poses[visited[i].first].residue(j+pair_1_helix_2_offset).atom_index("CA"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id1);
				
				core::id::AtomID const id2( poses[visited[i].second].residue(j+pair_2_helix_2_offset).atom_index("CA"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id2);
				atom_map[ id2 ] = id1;

				core::id::AtomID const id3( poses[visited[i].first].residue(j+pair_1_helix_2_offset).atom_index("C"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id3);
				core::id::AtomID const id4( poses[visited[i].second].residue(j+pair_2_helix_2_offset).atom_index("C"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id4);
				atom_map[ id4 ] = id3;

				core::id::AtomID const id5( poses[visited[i].first].residue(j+pair_1_helix_2_offset).atom_index("N"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id5);
				core::id::AtomID const id6( poses[visited[i].second].residue(j+pair_2_helix_2_offset).atom_index("N"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id6);
				atom_map[ id6 ] = id5;

				core::id::AtomID const id7( poses[visited[i].first].residue(j+pair_1_helix_2_offset).atom_index("O"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id7);
				core::id::AtomID const id8( poses[visited[i].second].residue(j+pair_2_helix_2_offset).atom_index("O"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id8);
				atom_map[ id8 ] = id7;
			}
			
			utility::vector1< numeric::xyzVector < core::Real > > helix_1_coords;
			utility::vector1< numeric::xyzVector < core::Real > > helix_2_coords;
			
			poses[visited[i].first].batch_get_xyz(helix_1_atoms, helix_1_coords);
			poses[visited[i].second].batch_get_xyz(helix_2_atoms, helix_2_coords);
			
			TR << "helix_1_coords size: " << helix_1_coords.size() << std::endl;
			TR << "helix_2_coords size: " << helix_2_coords.size() << std::endl;
			
			core::Real test_rmsd = numeric::model_quality::calc_rms(helix_1_coords,helix_2_coords);
			
			core::scoring::superimpose_pose(poses[visited[i].second], poses[visited[i].first]/*const*/, atom_map);

			TR << "RMSD for superimposition of pairs " << visited[i].first.print() << " " << visited[i].second.print() << ": " << test_rmsd << std::endl;

			if(master_pose.total_residue() == 0){
				master_pose = poses[visited[i].first];
			}
		}

		std::string frag_name = utility::to_string(tree_id) + "_" + utility::to_string(i) + "_" + visited[i].second.print() + ".pdb";
		output_poses.insert(std::make_pair(frag_name, poses[visited[i].second]));
		
		for(std::map<core::Size,core::Size>::const_iterator it=visited[i].second.helix_1_positions_.begin();
			it!=visited[i].second.helix_1_positions_.end(); ++it){
			if(visited[i].second.reversed()){
				TR << "ROUND: " << i << "-helix 1 mapping " << it->first+helix_size << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses[visited[i].second].residue(it->first+helix_size).name1());
			}
			else{
				TR << "ROUND: " << i << "-helix 1 mapping " << it->first << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses[visited[i].second].residue(it->first).name1());
			}
		}
		for(std::map<core::Size,core::Size>::const_iterator it=visited[i].second.helix_2_positions_.begin();
			it!=visited[i].second.helix_2_positions_.end(); ++it){
			if(visited[i].second.reversed()){
				TR << "ROUND: " << i << "-helix 2 mapping " << it->first << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses[visited[i].second].residue(it->first).name1());
			}
			else{
				TR << "ROUND: " << i << "-helix 2 mapping " << it->first+helix_size << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses[visited[i].second].residue(it->first+helix_size).name1());
			}
		}
	}
	master_pose.dump_pdb(output_name + ".pdb");
	
	dumpResfile(sequence_map, output_name + ".res");
	for(std::map<std::string, core::pose::Pose>::const_iterator out_it=output_poses.begin();
		out_it!= output_poses.end();
		++out_it){
		out_it->second.dump_pdb(out_it->first);
	}
}
