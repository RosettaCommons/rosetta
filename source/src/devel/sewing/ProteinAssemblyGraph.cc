// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProteinAssemblyGraph.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <devel/sewing/ProteinAssemblyGraph.hh>
#include <devel/sewing/util.hh>

//Core headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/Etable.hh>

//Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

//Basic headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

//Protocol headers
#include <protocols/features/ProteinSilentReport.hh>

//Boost


namespace devel {
namespace sewing {

static THREAD_LOCAL basic::Tracer TR( "ProteinAssesmblyGraph" );

ProteinAssemblyGraph::ProteinAssemblyGraph(
	utility::sql_database::sessionOP db_session
):
	db_session_(db_session),
	inter_structure_edge_counter_(0)
{
	scorefxn_ = core::scoring::get_score_function();
}

void
ProteinAssemblyGraph::addIntraStructureEdge(
	NodeOP node1,
	NodeOP node2
){
	nodes_.insert(node1);
	nodes_.insert(node2);
	AdjacencyMap::iterator adjacent = intra_structure_adjacency_map_.find(node1);

	if ( adjacent == intra_structure_adjacency_map_.end() ) {
		std::set<NodeOP> new_adjacency_list;
		intra_structure_adjacency_map_.insert(std::pair< NodeOP, std::set<NodeOP> >(node1,new_adjacency_list));
	}
	intra_structure_adjacency_map_[node1].insert(node2);
}

void
ProteinAssemblyGraph::addInterStructureEdge(
	NodeOP node1,
	NodeOP node2
){
	nodes_.insert(node1);
	nodes_.insert(node2);
	AdjacencyMap::iterator adjacent = inter_structure_adjacency_map_.find(node1);

	if ( adjacent == inter_structure_adjacency_map_.end() ) {
		std::set<NodeOP> new_adjacency_list;
		inter_structure_adjacency_map_.insert(std::pair< NodeOP, std::set<NodeOP> >(node1,new_adjacency_list));
	}
	inter_structure_adjacency_map_[node1].insert(node2);
	inter_structure_edge_counter_++;
}

utility::vector1<ProteinAssemblyGraph::EdgeList>
ProteinAssemblyGraph::find_all_trees(
	core::Size num_helices
){
	TR << "Total nodes in graph: " << nodes_.size() << std::endl;
	TR << "Total inter-bundle edges: " << total_inter_structure_edges() << std::endl;
	core::Real round=0.0;
	core::Real print_cutoff = 0.10;
	for ( std::set<NodeOP>::iterator it=nodes_.begin(); it!=nodes_.end(); ++it ) {
		utility::vector1<NodeOP> visited;
		visited.push_back((*it));

		//  DEBUG
		//  if((*it)->pair_id()==1052)
		//  {

		core::Size path_size = get_path_size(num_helices);

		ProteinAssemblyGraph::EdgeList visted_edges;
		ProteinAssemblyGraph::EdgeSet possible_inter_edges;
		ProteinAssemblyGraph::EdgeSet possible_intra_edges;
		tree_finder(visted_edges, visited, path_size, 0, false, possible_inter_edges, possible_intra_edges);

		round+=1.0;
		core::Real progress = round/nodes_.size();
		if ( progress > print_cutoff ) {
			TR << print_cutoff*100 << "% complete - " << tree_list_.size() << " total trees" << std::endl;
			print_cutoff+=0.10;
		}
		//  }
	}
	return tree_list_;
}


void
ProteinAssemblyGraph::print_tree_simple(
	EdgeList visited
){
	for ( core::Size i=1; i<=visited.size(); ++i ) {
		TR << "(" << visited[i].first->print() << ", " << visited[i].second->print() << ") ";
	}
	TR << std::endl;
}

core::Real
ProteinAssemblyGraph::bb_score(
	core::pose::Pose & pose,
	core::Size new_residues_begin,
	core::Size new_residues_end
){

	// score the bb-bb energy between chains
	// This part written by P.Doug Renfrew
	// make vectors of bb atoms in each chain individually
	// the master pose will always be chain 1.
	// need to make a vector of all atoms in the chain you are comparing too

	using namespace core;


	utility::vector1<core::conformation::Atom> new_helix_bb_atoms;
	utility::vector1<core::conformation::Atom> other_bb_atoms;

	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		core::conformation::Residue const & res( pose.residue(j) );
		core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
		for ( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
			if ( j>=new_residues_begin && j<=new_residues_end ) {
				new_helix_bb_atoms.push_back( res.atom(jj) );
			} else {
				other_bb_atoms.push_back( res.atom(jj) );
			}
		}
	}

	// get instance of etable energy method
	core::scoring::methods::EnergyMethodOptions const & emo(scorefxn_->energy_method_options());
	core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo).lock()));
	core::scoring::etable::AnalyticEtableEvaluator ete(et);

	// iterate over both sets of atom and add into one emapvector
	core::scoring::EMapVector tbemv;
	core::Real atr_wt( scorefxn_->get_weight(core::scoring::fa_atr) );
	core::Real rep_wt( scorefxn_->get_weight(core::scoring::fa_rep) );
	for ( Size ii = 1; ii <= new_helix_bb_atoms.size(); ++ii ) {
		for ( Size jj = 1; jj <= other_bb_atoms.size(); ++jj ) {
			//calc distr squared
			Real d2(new_helix_bb_atoms[ii].xyz().distance_squared(other_bb_atoms[jj].xyz()));
			ete.atom_pair_energy( new_helix_bb_atoms[ii], other_bb_atoms[jj], 1, tbemv, d2 );
		}
	}
	core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );

	return bb_energy;

}//end bb_score

void
ProteinAssemblyGraph::dump_resfile(
	std::map<core::Size, utility::vector1<char> > sequence_map,
	std::string filename
){

	utility::io::ozstream resfile;
	resfile.open(filename);

	resfile << "NATRO\n"
		"EX 1 EX 2\n"
		"USE_INPUT_SC\n"
		"start\n";
	for ( std::map<core::Size, utility::vector1<char> >::const_iterator it=sequence_map.begin();
			it != sequence_map.end(); ++it ) {
		core::Size position = it->first;
		utility::vector1<char> aas = it->second;

		resfile << position << " A PIKAA ";
		for ( core::Size i=1; i<=aas.size(); ++i ) {
			resfile << aas[i];
		}
		resfile << "\n";
	}
	resfile.close();
}

/// @brief Create poses for the given EdgeList.
void
ProteinAssemblyGraph::print_tree(
	EdgeList visited,
	core::Size tree_id
){

	std::string output_name = "tree_" + utility::to_string(tree_id) + "_";
	core::pose::Pose master_pose;
	utility::vector1<core::pose::Pose> master_poses;
	//core::Size repeat_counter=1;

	Size starttime = time( NULL );

	std::string select_pair_residues =
		"SELECT r.seqpos AS resnum\n"
		"FROM protein_residue_conformation r\n"
		"WHERE\n"
		"\tr.struct_id = ? AND"
		"\t(r.seqpos BETWEEN ? AND ? OR\n"
		"\tr.seqpos BETWEEN ? AND ?)\n";

	cppdb::statement select_pair_residues_stmt =
		basic::database::safely_prepare_statement(select_pair_residues, db_session_);

	protocols::features::ProteinSilentReportOP protein_silent_report( new protocols::features::ProteinSilentReport() );

	std::map<std::string, core::pose::Pose> output_poses;
	//std::map<Node, core::pose::Pose> poses;
	//retrieve all the structural info we need for each bundle pair in the given tree
	for ( core::Size i=1; i<=visited.size(); ++i ) {

		NodeOP cur_node;
		for ( core::Size j=1; j<=2; ++j ) {
			if ( j==1 ) {
				cur_node=visited[i].first;
			} else {
				cur_node=visited[i].second;
			}

			if ( poses_.find(cur_node) == poses_.end() ) { //new node

				std::set<core::Size> resnums;
				core::pose::Pose pose;

				select_pair_residues_stmt.bind(1, cur_node->struct_id());
				select_pair_residues_stmt.bind(2, cur_node->helix_1_begin());
				select_pair_residues_stmt.bind(3, cur_node->helix_1_end());
				select_pair_residues_stmt.bind(4, cur_node->helix_2_begin());
				select_pair_residues_stmt.bind(5, cur_node->helix_2_end());
				cppdb::result res = basic::database::safely_read_from_database(select_pair_residues_stmt);
				while ( res.next() ) {
					core::Size resnum;
					res >> resnum;
					resnums.insert(resnum);
				}

				protein_silent_report->load_pose(db_session_, cur_node->struct_id(), resnums, pose);
				poses_[cur_node]=pose;

			}
		}
	}
	Size endtime = time( NULL );
	TR << "Populated node poses in time: " << endtime - starttime << std::endl;

	//Clear tree memory from previous assemblies
	for ( core::Size i=1; i<=visited.size(); ++i ) {
		visited[i].first->helix_1_positions_.clear();
		visited[i].second->helix_1_positions_.clear();
		visited[i].first->helix_2_positions_.clear();
		visited[i].second->helix_2_positions_.clear();
	}

	//Assemble into the final bundle
	std::map<core::Size, utility::vector1<char> > sequence_map;//for resfile
	NativeRotamersMap native_residue_map;//for resfile
	std::set<NodeOP> printed_nodes;
	bool repeat=false;
	for ( core::Size i=1; i<=visited.size(); ++i ) {

		//For repeat proteins, clear the printed nodes list at the end of each repeat
		if ( visited[i].second == visited[1].first ) {
			printed_nodes.clear();
		}

		//temp for 4 helix bundles only
		if ( visited.size() >= i+1 && visited[i+1].second == visited[1].first ) {
			repeat=true;
		}

		NodeOP edge_start = visited[i].first;
		NodeOP edge_end = visited[i].second;
		core::Size helix_size = edge_start->helix_1_end()-edge_start->helix_1_begin()+1;

		//For the first edge dump the starting node to a pdb and then map each residue in the node to
		//to a residue number in the assembled pdb (this is the helix_1/2_positions map)
		if ( i==1 ) {
			std::string frag_name = utility::to_string(tree_id) + "_" + utility::to_string(i) + "_" + edge_start->print() + ".pdb";
			output_poses.insert(std::make_pair(frag_name, poses_[edge_start]));

			//populate the residue mapping for the first two helices
			core::Size position_counter=0;
			if ( edge_start->helix_2_begin() > edge_start->helix_1_begin() ) {
				for ( core::Size j=1; j<=helix_size; ++j ) {
					position_counter++;
					edge_start->helix_1_positions_[j]=position_counter;
				}
				for ( core::Size j=1; j<=helix_size; ++j ) {
					position_counter++;
					edge_start->helix_2_positions_[j]=position_counter;
				}
			} else { //helices will be assembled in the opposite sequence order they appear in the native
				for ( core::Size j=1; j<=helix_size; ++j ) {
					position_counter++;
					edge_start->helix_2_positions_[j]=position_counter;
				}
				for ( core::Size j=1; j<=helix_size; ++j ) {
					position_counter++;
					edge_start->helix_1_positions_[j]=position_counter;
				}
			}

			//Go through the residue mapping and
			for ( std::map<core::Size,core::Size>::const_iterator it=edge_start->helix_1_positions_.begin();
					it!=edge_start->helix_1_positions_.end(); ++it ) {
				core::Size helix_resnum = it->first;
				core::Size assembled_pose_resnum = it->second;

				if ( edge_start->reversed() ) {
					sequence_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum+helix_size).name1());
					native_residue_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum+helix_size).clone());
				} else {
					sequence_map[it->second].push_back(poses_[edge_start].residue(helix_resnum).name1());
					native_residue_map[it->second].push_back(poses_[edge_start].residue(helix_resnum).clone());
				}
			}
			for ( std::map<core::Size,core::Size>::const_iterator it=edge_start->helix_2_positions_.begin();
					it!=edge_start->helix_2_positions_.end(); ++it ) {
				core::Size helix_resnum = it->first;
				core::Size assembled_pose_resnum = it->second;

				if ( edge_start->reversed() ) {
					sequence_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum).name1());
					native_residue_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum).clone());
				} else {
					sequence_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum+helix_size).name1());
					native_residue_map[assembled_pose_resnum].push_back(poses_[edge_start].residue(helix_resnum+helix_size).clone());
				}
			}
		}

		//Pairs come from the same bundle, only superimpose the matching helix
		if ( edge_start->bundle_id() == edge_end->bundle_id() ) {
			core::Size first_pair_offset;
			core::Size second_pair_offset;
			core::Size new_helix_start;
			//all helices should be the same size

			//Find which helix matches between the two pairs
			if ( edge_start->helix_1_id()== edge_end->helix_1_id() ) {
				first_pair_offset=(edge_start->helix_1_begin()< edge_start->helix_2_begin()) ? 0 : helix_size;
				second_pair_offset=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? 0 : helix_size;
				new_helix_start=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? helix_size+1 : 1;

				edge_end->helix_1_positions_= edge_start->helix_1_positions_;
				core::Size position_start = master_pose.total_residue();
				for ( core::Size j=1; j<=helix_size; ++j ) {
					edge_end->helix_2_positions_[j]=position_start+j;
				}
			} else if ( edge_start->helix_1_id()== edge_end->helix_2_id() ) {
				first_pair_offset=(edge_start->helix_1_begin()< edge_start->helix_2_begin()) ? 0 : helix_size;
				second_pair_offset=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? helix_size : 0;
				new_helix_start=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? 1 : helix_size+1;

				edge_end->helix_2_positions_= edge_start->helix_1_positions_;
				core::Size position_start = master_pose.total_residue();
				for ( core::Size j=1; j<=helix_size; ++j ) {
					edge_end->helix_1_positions_[j]=position_start+j;
				}
			} else if ( edge_start->helix_2_id()== edge_end->helix_1_id() ) {
				first_pair_offset=(edge_start->helix_1_begin()< edge_start->helix_2_begin()) ? helix_size : 0;
				second_pair_offset=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? 0 : helix_size;
				new_helix_start=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? helix_size+1 : 1;

				edge_end->helix_1_positions_= edge_start->helix_2_positions_;
				core::Size position_start = master_pose.total_residue();
				for ( core::Size j=1; j<=helix_size; ++j ) {
					edge_end->helix_2_positions_[j]=position_start+j;
				}
			} else if ( edge_start->helix_2_id()== edge_end->helix_2_id() ) {
				first_pair_offset=(edge_start->helix_1_begin()< edge_start->helix_2_begin()) ? helix_size : 0;
				second_pair_offset=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? helix_size : 0;;
				new_helix_start=(edge_end->helix_1_begin()< edge_end->helix_2_begin()) ? 1 : helix_size+1;

				edge_end->helix_2_positions_= edge_start->helix_2_positions_;
				core::Size position_start = master_pose.total_residue();
				for ( core::Size j=1; j<=helix_size; ++j ) {
					edge_end->helix_1_positions_[j]=position_start+j;
				}
			} else {
				utility_exit_with_message(
					"No helix ids shared between node " + utility::to_string(edge_start->pair_id())+" and "+
					utility::to_string(edge_end->pair_id())+". Something has gone wrong with helical assembly!");
				//Gets rid of some warnings
				new_helix_start=0;
				first_pair_offset=0;
				second_pair_offset=0;
			}

			core::id::AtomID_Map< core::id::AtomID > atom_map;
			atom_map.clear();
			core::pose::initialize_atomid_map( atom_map, poses_[edge_end], core::id::BOGUS_ATOM_ID );

			for ( core::Size j=1; j<= helix_size; j++ ) {
				core::id::AtomID const id1( poses_[edge_start].residue(j+first_pair_offset).atom_index("CA"), j+first_pair_offset);
				core::id::AtomID const id2( poses_[edge_end].residue(j+second_pair_offset).atom_index("CA"), j+second_pair_offset);
				atom_map[ id2 ] = id1;

				core::id::AtomID const id3( poses_[edge_start].residue(j+first_pair_offset).atom_index("C"), j+first_pair_offset);
				core::id::AtomID const id4( poses_[edge_end].residue(j+second_pair_offset).atom_index("C"), j+second_pair_offset);
				atom_map[ id4 ] = id3;

				core::id::AtomID const id5( poses_[edge_start].residue(j+first_pair_offset).atom_index("N"), j+first_pair_offset);
				core::id::AtomID const id6( poses_[edge_end].residue(j+second_pair_offset).atom_index("N"), j+second_pair_offset);
				atom_map[ id6 ] = id5;

				core::id::AtomID const id7( poses_[edge_start].residue(j+first_pair_offset).atom_index("O"), j+first_pair_offset);
				core::id::AtomID const id8( poses_[edge_end].residue(j+second_pair_offset).atom_index("O"), j+second_pair_offset);
				atom_map[ id8 ] = id7;
			}
			core::scoring::superimpose_pose(poses_[edge_end], poses_[edge_start]/*const*/, atom_map);

			//don't print the same node more than once
			if ( printed_nodes.find(edge_start) == printed_nodes.end() ) {
				if ( repeat ) {
					core::pose::append_subpose_to_pose(master_pose, poses_[edge_end], new_helix_start, new_helix_start+helix_size-1, true);
					repeat=false;
				} else {
					core::pose::append_subpose_to_pose(master_pose, poses_[edge_end], new_helix_start, new_helix_start+helix_size-1, false);
				}
				printed_nodes.insert(edge_start);
			}

			//   for(core::Size i=1; i<= master_poses.size(); ++i)
			//   {
			//    core::pose::append_subpose_to_pose(master_poses[i], poses_[edge_end], new_helix_start, new_helix_start+helix_size-1, false);
			//    core::pose::append_subpose_to_pose(master_poses[i], poses_[edge_end], new_helix_start, new_helix_start+helix_size-1, false);
			//   }
			output_name+=utility::to_string(edge_start->pair_id()) + "_";

			//   core::Real clash_score = bb_score(master_pose, pose_size_before_addition+1, master_pose.total_residue());
			//  TR << "Backbone-backbone score for addition of pair " << tag_1 << ": " << clash_score << std::endl;
			//   if(clash_score > 5){
			//    return;
			//   }

		} else {
			//Pairs come from different bundles, superimpose both helices
			edge_end->helix_1_positions_= edge_start->helix_1_positions_;
			edge_end->helix_2_positions_= edge_start->helix_2_positions_;

			core::id::AtomID_Map< core::id::AtomID > atom_map;
			atom_map.clear();
			core::pose::initialize_atomid_map( atom_map, poses_[edge_end], core::id::BOGUS_ATOM_ID );

			utility::vector1<core::id::AtomID> helix_1_atoms;
			utility::vector1<core::id::AtomID> helix_2_atoms;

			core::Size pair_1_helix_1_offset;
			core::Size pair_1_helix_2_offset;
			core::Size pair_2_helix_1_offset;
			core::Size pair_2_helix_2_offset;
			if ( edge_start->helix_1_begin() < edge_start->helix_2_begin() ) {
				pair_1_helix_1_offset = 0;
				pair_1_helix_2_offset = helix_size;
			} else {
				pair_1_helix_1_offset = helix_size;
				pair_1_helix_2_offset = 0;
			}

			if ( edge_end->helix_1_begin() < edge_end->helix_2_begin() ) {
				pair_2_helix_1_offset = 0;
				pair_2_helix_2_offset = helix_size;
			} else {
				pair_2_helix_1_offset = helix_size;
				pair_2_helix_2_offset = 0;
			}

			for ( core::Size j=1; j<=helix_size; j++ ) {
				core::id::AtomID const id1( poses_[edge_start].residue(j+pair_1_helix_1_offset).atom_index("CA"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id1);

				core::id::AtomID const id2( poses_[edge_end].residue(j+pair_2_helix_1_offset).atom_index("CA"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id2);
				atom_map[ id2 ] = id1;

				core::id::AtomID const id3( poses_[edge_start].residue(j+pair_1_helix_1_offset).atom_index("C"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id3);
				core::id::AtomID const id4( poses_[edge_end].residue(j+pair_2_helix_1_offset).atom_index("C"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id4);
				atom_map[ id4 ] = id3;

				core::id::AtomID const id5( poses_[edge_start].residue(j+pair_1_helix_1_offset).atom_index("N"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id5);
				core::id::AtomID const id6( poses_[edge_end].residue(j+pair_2_helix_1_offset).atom_index("N"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id6);
				atom_map[ id6 ] = id5;

				core::id::AtomID const id7( poses_[edge_start].residue(j+pair_1_helix_1_offset).atom_index("O"), j+pair_1_helix_1_offset);
				helix_1_atoms.push_back(id7);
				core::id::AtomID const id8( poses_[edge_end].residue(j+pair_2_helix_1_offset).atom_index("O"), j+pair_2_helix_1_offset);
				helix_2_atoms.push_back(id8);
				atom_map[ id8 ] = id7;
			}

			for ( core::Size j=1; j<=helix_size; j++ ) {
				core::id::AtomID const id1( poses_[edge_start].residue(j+pair_1_helix_2_offset).atom_index("CA"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id1);

				core::id::AtomID const id2( poses_[edge_end].residue(j+pair_2_helix_2_offset).atom_index("CA"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id2);
				atom_map[ id2 ] = id1;

				core::id::AtomID const id3( poses_[edge_start].residue(j+pair_1_helix_2_offset).atom_index("C"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id3);
				core::id::AtomID const id4( poses_[edge_end].residue(j+pair_2_helix_2_offset).atom_index("C"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id4);
				atom_map[ id4 ] = id3;

				core::id::AtomID const id5( poses_[edge_start].residue(j+pair_1_helix_2_offset).atom_index("N"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id5);
				core::id::AtomID const id6( poses_[edge_end].residue(j+pair_2_helix_2_offset).atom_index("N"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id6);
				atom_map[ id6 ] = id5;

				core::id::AtomID const id7( poses_[edge_start].residue(j+pair_1_helix_2_offset).atom_index("O"), j+pair_1_helix_2_offset);
				helix_1_atoms.push_back(id7);
				core::id::AtomID const id8( poses_[edge_end].residue(j+pair_2_helix_2_offset).atom_index("O"), j+pair_2_helix_2_offset);
				helix_2_atoms.push_back(id8);
				atom_map[ id8 ] = id7;
			}

			utility::vector1< numeric::xyzVector < core::Real > > helix_1_coords;
			utility::vector1< numeric::xyzVector < core::Real > > helix_2_coords;

			poses_[edge_start].batch_get_xyz(helix_1_atoms, helix_1_coords);
			poses_[edge_end].batch_get_xyz(helix_2_atoms, helix_2_coords);

			TR << "helix_1_coords size: " << helix_1_coords.size() << std::endl;
			TR << "helix_2_coords size: " << helix_2_coords.size() << std::endl;

			core::Real test_rmsd = numeric::model_quality::calc_rms(helix_1_coords,helix_2_coords);

			core::scoring::superimpose_pose(poses_[edge_end], poses_[edge_start]/*const*/, atom_map);

			TR << "RMSD for superimposition of pairs " << edge_start->print() << " " << edge_end->print() << ": " << test_rmsd << std::endl;

			if ( master_pose.total_residue() == 0 ) {
				master_pose=print_first(visited, poses_);
				//    master_poses.push_back(poses_[edge_start])
				//    master_poses.push_back(poses_[edge_end])
			} else {
				//need to replace the two matching helices with new ones in all the master_poses
			}
		}

		std::string frag_name = utility::to_string(tree_id) + "_" + utility::to_string(i) + "_" + edge_end->print() + ".pdb";
		output_poses.insert(std::make_pair(frag_name, poses_[edge_end]));

		for ( std::map<core::Size,core::Size>::const_iterator it=edge_end->helix_1_positions_.begin();
				it!=edge_end->helix_1_positions_.end(); ++it ) {
			if ( edge_end->reversed() ) {
				//    TR << "ROUND: " << i << "-helix 1 mapping " << it->first+helix_size << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses_[edge_end].residue(it->first+helix_size).name1());
				native_residue_map[it->second].push_back(poses_[edge_end].residue(it->first+helix_size).clone());
			} else {
				//    TR << "ROUND: " << i << "-helix 1 mapping " << it->first << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses_[edge_end].residue(it->first).name1());
				native_residue_map[it->second].push_back(poses_[edge_end].residue(it->first).clone());
			}
		}
		for ( std::map<core::Size,core::Size>::const_iterator it=edge_end->helix_2_positions_.begin();
				it!=edge_end->helix_2_positions_.end(); ++it ) {
			if ( edge_end->reversed() ) {
				//    TR << "ROUND: " << i << "-helix 2 mapping " << it->first << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses_[edge_end].residue(it->first).name1());
				native_residue_map[it->second].push_back(poses_[edge_end].residue(it->first).clone());
			} else {
				//    TR << "ROUND: " << i << "-helix 2 mapping " << it->first+helix_size << " to bundle position " << it->second << std::endl;
				sequence_map[it->second].push_back(poses_[edge_end].residue(it->first+helix_size).name1());
				native_residue_map[it->second].push_back(poses_[edge_end].residue(it->first+helix_size).clone());
			}
		}
	}
	// master_pose.dump_pdb(output_name + ".pdb");
	master_pose.dump_pdb("tree_" + utility::to_string(tree_id) + ".pdb");
	//std::string repeat_name = "tree_" + utility::to_string(tree_id) + "_repeat_" + utility::to_string(repeat_counter);

	dump_resfile(sequence_map, output_name + ".res");
	devel::sewing::dump_native_residue_file(native_residue_map, output_name + ".rots");
	for ( std::map<std::string, core::pose::Pose>::const_iterator out_it=output_poses.begin();
			out_it!= output_poses.end();
			++out_it ) {
		out_it->second.dump_pdb(out_it->first);
	}
}

core::pose::Pose
ProteinAssemblyGraph::print_first(
	EdgeList const & visited_edges,
	std::map<NodeOP, core::pose::Pose> const & poses
){
	return poses.find(visited_edges[1].first)->second;
}

std::set<NodeOP>
ProteinAssemblyGraph::nodes()
{
	return nodes_;
}

core::Size
ProteinAssemblyGraph::total_inter_structure_edges()
{
	return inter_structure_edge_counter_;
}

utility::vector1<ProteinAssemblyGraph::EdgeList>
ProteinAssemblyGraph::tree_list()
{
	return tree_list_;
}

} //sewing namespace
} //devel namespace
