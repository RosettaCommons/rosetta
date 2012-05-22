// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixAssemblyGraph.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <devel/helixAssembly/HelixAssemblyGraph.hh>

//Core headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
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

//Basic headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

//Protocol headers
#include <protocols/features/ProteinSilentReport.hh>

static basic::Tracer TR("HelixAssesmblyGraph");

//****Node class****//

Node::Node(core::Size pair_id):
pair_id_(pair_id)
{}

core::Size Node::pair_id() const{
	return pair_id_;
}

std::string Node::print() const{
	return utility::to_string(pair_id_);
}

bool Node::operator<(const Node & rhs) const{
	return (pair_id_ < rhs.pair_id());
}

bool Node::operator==(const Node & rhs) const{
	return (pair_id_ == rhs.pair_id());
}

//****Graph class****//

HelixAssemblyGraph::HelixAssemblyGraph(utility::sql_database::sessionOP db_session, bool simple_print):
db_session_(db_session),
simple_print_(simple_print),
inter_bundle_edge_counter_(0)
{}

void HelixAssemblyGraph::addIntraBundleEdge(Node node1, Node node2) {
	nodes_.insert(node1);
	nodes_.insert(node2);
	AdjacencyMap::iterator adjacent = intra_bundle_adjacency_map_.find(node1);
	
	if(adjacent == intra_bundle_adjacency_map_.end()){
		utility::vector1<Node> new_adjacency_list;
		intra_bundle_adjacency_map_.insert(std::pair< Node, utility::vector1<Node> >(node1,new_adjacency_list));
	}
	intra_bundle_adjacency_map_[node1].push_back(node2);
}

void HelixAssemblyGraph::addInterBundleEdge(Node node1, Node node2) {
	nodes_.insert(node1);
	nodes_.insert(node2);
	AdjacencyMap::iterator adjacent = inter_bundle_adjacency_map_.find(node1);
	
	if(adjacent == inter_bundle_adjacency_map_.end()){
		utility::vector1<Node> new_adjacency_list;
		inter_bundle_adjacency_map_.insert(std::pair< Node, utility::vector1<Node> >(node1,new_adjacency_list));
	}
	inter_bundle_adjacency_map_[node1].push_back(node2);
	inter_bundle_edge_counter_++;
}

//Depth first search prints all paths of a given size
void HelixAssemblyGraph::depthFirst(utility::vector1<Node> visited, size_t path_size){

	TR.Debug << "Dept first called. visited size is: " << visited.size() << std::endl;

	//check to see if we should use inter or intra bundle edges and find adjacent nodes accordingly
	utility::vector1<Node> nodes;
	if(visited.size() % 2 == 0){
		nodes = intra_bundle_adjacency_map_[visited[visited.size()]];
	}
	else{
		nodes = inter_bundle_adjacency_map_[visited[visited.size()]];
	}
	
	// examine adjacent nodes for end condition (path length in my case)
	for(core::Size i=1; i<=nodes.size(); ++i){
		if (visited.contains(nodes[i])){
			continue;
		}

		//Check for end condition (path size in this case. could also be cycle, etc.)
		if (visited.size()+1 == path_size) {
			visited.push_back(nodes[i]);
			if(simple_print_){
				printPathSimple(visited);
			}
			else{
				printPath(visited);
			}
			path_list_.push_back(visited);
			visited.pop_back();
			break;
		}
	}
	
	for(core::Size i=1; i<=nodes.size(); ++i){
		if (visited.contains(nodes[i])) {
			continue;
		}
		visited.push_back(nodes[i]);
		depthFirst(visited, path_size);
		visited.pop_back();
	}
}

void HelixAssemblyGraph::printPathSimple(utility::vector1<Node> visited){
	for(core::Size i=1; i<=visited.size(); ++i){
		TR << visited[i].print() << " ";
	}
	TR << std::endl;
}

void HelixAssemblyGraph::printPath(utility::vector1<Node> visited){
	
	std::string output_name="bundle_" + utility::to_string(path_list_.size()) + "_";
	core::pose::Pose master_pose;
	
	std::string select_pairs = 
	"SELECT s.struct_id AS struct_id, r.seqpos AS resnum, bp.pair_id AS tag,\n"
	"	bp.helix_id_1, bp.helix_1_begin, bp.helix_1_end,\n"
	"	bp.helix_id_2, bp.helix_2_begin, bp.helix_2_end\n"
	"FROM bundle_pairs bp\n"
	"JOIN helix_bundles hb ON\n"
	"	bp.bundle_id=hb.bundle_id\n"
	"JOIN structures s ON\n"
	"	hb.struct_id = s.struct_id\n"
	"JOIN protein_residue_conformation r ON\n"
	"	(r.seqpos BETWEEN bp.helix_1_begin AND bp.helix_1_end OR\n"
	"	r.seqpos BETWEEN bp.helix_2_begin AND bp.helix_2_end) AND\n"
	"	r.struct_id = s.struct_id\n"
	"WHERE bp.pair_id=?";
	
	cppdb::statement select_pairs_stmt = basic::database::safely_prepare_statement(select_pairs, db_session_);
	
	protocols::features::ProteinSilentReportOP protein_silent_report = new protocols::features::ProteinSilentReport();
	
	core::pose::Pose prev_pose;
	core::Size prev_helix_id_1;	
	core::Size prev_helix_1_begin;
	core::Size prev_helix_1_end;
	
	core::Size prev_helix_id_2;
	core::Size prev_helix_2_begin;
	core::Size prev_helix_2_end;
	
	//go two by two aligning the bundles that match together
	for(core::Size i=1; i<=visited.size(); i+=2){
		
		core::pose::Pose pose_1;
		boost::uuids::uuid struct_id_1;
		std::string tag_1;
		std::set<core::Size> resnums_1;
		core::Size pair_1_helix_id_1;
		core::Size pair_1_helix_1_begin;
		core::Size pair_1_helix_1_end;
		
		core::Size pair_1_helix_id_2;
		core::Size pair_1_helix_2_begin;
		core::Size pair_1_helix_2_end;
		
		select_pairs_stmt.bind(1, visited[i].pair_id());
		cppdb::result res = basic::database::safely_read_from_database(select_pairs_stmt);
		while(res.next()){
			core::Size resnum;
			res >> struct_id_1 >> resnum >> tag_1 >>
			pair_1_helix_id_1 >> pair_1_helix_1_begin >> pair_1_helix_1_end >>
			pair_1_helix_id_2 >> pair_1_helix_2_begin >> pair_1_helix_2_end;
			resnums_1.insert(resnum);
		}
		protein_silent_report->load_pose(db_session_, struct_id_1, resnums_1, pose_1);
		
		//superimpose matching helix from pose_1 onto prev_pose
		if(prev_pose.total_residue()>0){
			core::Size pair_1_offset;
			core::Size prev_pair_offset;
			core::Size new_helix_start;
			//all helices should be the same size
			core::Size helix_size = pair_1_helix_1_end-pair_1_helix_1_begin+1;
			if(pair_1_helix_id_1 == prev_helix_id_1){
				pair_1_offset=0;
				prev_pair_offset=0;
				new_helix_start=helix_size;
			}
			else if(pair_1_helix_id_1 == prev_helix_id_2){
				pair_1_offset=0;
				prev_pair_offset=helix_size;
				new_helix_start=helix_size;
			}
			else if(pair_1_helix_id_2 == prev_helix_id_1){
				pair_1_offset=helix_size;
				prev_pair_offset=0;
				new_helix_start=1;
			}
			else if(pair_1_helix_id_2 == prev_helix_id_2){
				pair_1_offset=helix_size;
				prev_pair_offset=helix_size;
				new_helix_start=1;
			}
			else{
				utility_exit_with_message("No helix ids matched. something has gone wrong with helical assembly!");
			}
			
			core::id::AtomID_Map< core::id::AtomID > atom_map;
			atom_map.clear();
			core::pose::initialize_atomid_map( atom_map, pose_1, core::id::BOGUS_ATOM_ID );
			
			for(core::Size j=1; j<= helix_size; j++){
				core::id::AtomID const id1( prev_pose.residue(j+prev_pair_offset).atom_index("CA"), j+prev_pair_offset);
				core::id::AtomID const id2( pose_1.residue(j+pair_1_offset).atom_index("CA"), j+pair_1_offset);
				atom_map[ id2 ] = id1;
				
				core::id::AtomID const id3( prev_pose.residue(j+prev_pair_offset).atom_index("C"), j+prev_pair_offset);
				core::id::AtomID const id4( pose_1.residue(j+pair_1_offset).atom_index("C"), j+pair_1_offset);
				atom_map[ id4 ] = id3;
				
				core::id::AtomID const id5( prev_pose.residue(j+prev_pair_offset).atom_index("N"), j+prev_pair_offset);
				core::id::AtomID const id6( pose_1.residue(j+pair_1_offset).atom_index("N"), j+pair_1_offset);
				atom_map[ id6 ] = id5;
				
				core::id::AtomID const id7( prev_pose.residue(j+prev_pair_offset).atom_index("O"), j+prev_pair_offset);
				core::id::AtomID const id8( pose_1.residue(j+pair_1_offset).atom_index("O"), j+pair_1_offset);
				atom_map[ id8 ] = id7;
			}
			core::scoring::superimpose_pose(pose_1, prev_pose/*const*/, atom_map);
			
			combinePoses(master_pose, pose_1);
//			if(bb_score(master_pose, new_helix_start, new_helix_start+helix_size-1) > 5){
//				//BACKBONE CLASH FAIL
//				return;
//			}
		}
		
		core::pose::Pose pose_2;
		boost::uuids::uuid struct_id_2;
		std::string tag_2;
		std::set<core::Size> resnums_2;
		core::Size pair_2_helix_id_1;
		core::Size pair_2_helix_1_begin;
		core::Size pair_2_helix_1_end;
		
		core::Size pair_2_helix_id_2;
		core::Size pair_2_helix_2_begin;
		core::Size pair_2_helix_2_end;
		
		select_pairs_stmt.bind(1, visited[i+1].pair_id());
		res = basic::database::safely_read_from_database(select_pairs_stmt);
		while(res.next()){
			core::Size resnum;
			res >> struct_id_2 >> resnum >> tag_2 >>
			pair_2_helix_id_1 >> pair_2_helix_1_begin >> pair_2_helix_1_end >>
			pair_2_helix_id_2 >> pair_2_helix_2_begin >> pair_2_helix_2_end;
			resnums_2.insert(resnum);
		}
		protein_silent_report->load_pose(db_session_, struct_id_2, resnums_2, pose_2);
		
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		atom_map.clear();
		core::pose::initialize_atomid_map( atom_map, pose_2, core::id::BOGUS_ATOM_ID );
		
		for(core::Size j=1; j<= pose_2.total_residue(); j++){
			core::id::AtomID const id1( pose_1.residue(j).atom_index("CA"), j);
			core::id::AtomID const id2( pose_2.residue(j).atom_index("CA"), j );
			atom_map[ id2 ] = id1;
			
			core::id::AtomID const id3( pose_1.residue(j).atom_index("C"), j);
			core::id::AtomID const id4( pose_2.residue(j).atom_index("C"), j);
			atom_map[ id4 ] = id3;
			
			core::id::AtomID const id5( pose_1.residue(j).atom_index("N"), j);
			core::id::AtomID const id6( pose_2.residue(j).atom_index("N"), j);
			atom_map[ id6 ] = id5;
			
			core::id::AtomID const id7( pose_1.residue(j).atom_index("O"), j);
			core::id::AtomID const id8( pose_2.residue(j).atom_index("O"), j);
			atom_map[ id8 ] = id7;
		}
		core::scoring::superimpose_pose(pose_2, pose_1/*const*/, atom_map);
		pose_1.dump_pdb("bundle_" + utility::to_string(path_list_.size()) + "_pair_" + tag_1 + ".pdb");
		pose_2.dump_pdb("bundle_" + utility::to_string(path_list_.size()) + "_pair_" + tag_2 + ".pdb");
		
		combinePoses(master_pose, pose_2);
		
		output_name+=tag_1+"_"+tag_2+"_";
		
		prev_pose=pose_2;
		prev_helix_id_1=pair_2_helix_id_1;
		prev_helix_1_begin=pair_2_helix_1_begin;
		prev_helix_1_end=pair_2_helix_1_end;
		
		prev_helix_id_2=pair_2_helix_id_2;
		prev_helix_2_begin=pair_2_helix_2_begin;
		prev_helix_2_end=pair_2_helix_2_end;
	}
	output_name+=".pdb";
	master_pose.dump_pdb(output_name);
}

void HelixAssemblyGraph::combinePoses(core::pose::Pose & pose1, core::pose::Pose const & pose2){
	if(pose2.total_residue()>=1){
		pose1.append_residue_by_jump(pose2.residue(1), pose1.total_residue() , "", "", true/*start new chain*/);
		for(core::Size i=2; i<=pose2.total_residue(); ++i){
			if(pose2.residue(i).is_lower_terminus()){
				pose1.append_residue_by_jump(pose2.residue(i), pose1.total_residue(), "","", true);
			}
			else{
				pose1.append_residue_by_bond(pose2.residue(i));
			}
		}
	}
}

core::Real HelixAssemblyGraph::bb_score(core::pose::Pose & pose, core::Size new_helix_begin, core::Size new_helix_end){
	
	// score the bb-bb energy between chains
	// This part written by P.Doug Renfrew
	// make vectors of bb atoms in each chain individually
	// the master pose will always be chain 1.
	// need to make a vector of all atoms in the chain you are comparing too
	
	using namespace core;
	
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	
	utility::vector1<core::conformation::Atom> new_helix_bb_atoms;
	utility::vector1<core::conformation::Atom> other_bb_atoms;
	
	for( Size j = 1; j <= pose.total_residue(); ++j ) {
		core::conformation::Residue const & res( pose.residue(j) );
		core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
		for( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
			if( j>=new_helix_begin && j<=new_helix_end ){
				new_helix_bb_atoms.push_back( res.atom(jj) );
			}
			
			else{
				other_bb_atoms.push_back( res.atom(jj) );
			}
		}
	}
	
	//NOW SCORE!
	// get instance of etable energy method
	core::scoring::methods::EnergyMethodOptions const & emo(scorefxn->energy_method_options());
	core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo.etable_type())));
	core::scoring::etable::EtableEnergy ete( et, emo );
	
	// iterate over both sets of atom and add into one emapvector
	core::scoring::EMapVector tbemv;
	core::Real atr_wt( scorefxn->get_weight(core::scoring::fa_atr) );
	core::Real rep_wt( scorefxn->get_weight(core::scoring::fa_rep) );
	for ( Size ii = 1; ii <= new_helix_bb_atoms.size(); ++ii ){
		for ( Size jj = 1; jj <= other_bb_atoms.size(); ++jj ) {
			//calc distr squared
			Real d2(new_helix_bb_atoms[ii].xyz().distance_squared(other_bb_atoms[jj].xyz()));
			ete.atom_pair_energy( new_helix_bb_atoms[ii], other_bb_atoms[jj], 1, tbemv, d2 );
		}
	}
	core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );
	
	TR << "Backbone-backbone score: " << bb_energy << std::endl;
	return bb_energy;
	//return bb_energy;
}//end bb_score


std::set<Node> HelixAssemblyGraph::nodes(){
	return nodes_;
}

core::Size HelixAssemblyGraph::total_inter_bundle_edges(){
	return inter_bundle_edge_counter_;
}

//HelixAssemblyGraph::create_bundles(){
//	
//			
//	using namespace boost;
//	typedef adjacency_list<
//	vecS, //Use a vector for the vertex list
//	setS, //Use a set for the edge list (enforces no multi-graphs)
//	undirectedS, //This graph is undirected (structure similarity is undirected)
//	property<vertex_name_t, std::string>,
//	property<edge_weight_t, core::Real>//properties
//	> Graph;
//	
//	Graph g;
//	
//	property_map<Graph, vertex_name_t>::type vertex_id = get(vertex_name, g);
//	property_map<Graph, edge_weight_t>::type edge_rmsd = get(edge_weight, g);
//
//	typedef graph_traits<Graph>::vertex_descriptor Vertex;
//	typedef std::map<core::Size, Vertex> PairIdVertexMap;
//	PairIdVertexMap bundle_pairs;
//
////while(rex.next()){
//
//	PairIdVertexMap::iterator pos; 
//	bool inserted;
//	Vertex u, v;
//	
//	tie(pos, inserted) = bundle_pairs.insert(std::make_pair(1/*pair_id_1*/, Vertex()));
//	if (inserted) {//only add a vertex if we haven't added a vertex for this bundle pair
//		u = add_vertex(g); //add the new vertex
//		vertex_id[u] = 1/*vertex_id_1*/;
//		pos->second = u; //update the bundle pairs map to point to the newly added vertex
//	} else{
//		u = pos->second; //set the vertex to the one already in the graph
//	}
//	
//	tie(pos, inserted) = bundle_pairs.insert(std::make_pair(1/*pair_id_2*/, Vertex()));
//	if (inserted) {//only add a vertex if we haven't added a vertex for this bundle pair
//		v = add_vertex(g);
//		vertex_id[v] = 2/*vertex_id_2*/;
//		pos->second = v;
//	} else{
//		v = pos->second;
//	}
//	
//	graph_traits<Graph>::edge_descriptor e;
//	tie(e, inserted) = add_edge(u, v, g);
//	if (inserted){
//		edge_rmsd[e] = 5.3/*rmsd*/;
//	}
////}
	

	//Select all the vertices (bundle pairs) and edges (<0.5rmsd bundle_comparisons)	
//}