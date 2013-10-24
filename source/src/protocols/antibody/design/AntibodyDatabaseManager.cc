// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDatabaseManager.cc
/// @brief Handles all loading of CDR, Framework, and cluster/dmap info from external database file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Antibody Headers
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/CDRClusterEnumManager.hh>


//Core Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AA.hh>

//Database Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <protocols/features/ProteinSilentReport.hh>
#include <basic/database/sql_utils.hh>
#include <boost/uuid/uuid.hpp>
#include <cppdb/frontend.h>
#include <basic/database/open.hh>

//Option Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/inout.OptionKeys.gen.hh>

//Basic Headers
#include <basic/Tracer.hh>
#include <string>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/PyAssert.hh>
#include <boost/algorithm/string.hpp>
#include <fstream>


static basic::Tracer TR("antibody.design.AntibodyDatabaseManager");

namespace protocols {
namespace antibody {
namespace design {

using std::string;
using std::endl;
using core::pose::Pose;
using utility::file::FileName;
using utility::vector1;
using utility::sql_database::DatabaseSessionManager;
using protocols::features::ProteinSilentReport;
using namespace protocols::features;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::antibody;

typedef std::map< CDRNameEnum, CDRGraftInstructions > GraftInstructions;
typedef std::map< CDRNameEnum, CDRDesignInstructions > DesignInstructions;
typedef std::map< CDRNameEnum, vector1<PoseOP> > CDRSet;
typedef std::map< CDRNameEnum, vector1< CDRClusterEnum > > CDRClusterMap;
typedef std::map< CDRNameEnum, vector1< std::string > > PDBMap;

AntibodyDatabaseManager::AntibodyDatabaseManager(){
	//protein_silent_report_ = new protocols::features::ProteinSilentReport();
	string default_path = basic::database::full_name(basic::options::option[basic::options::OptionKeys::antibody::design::antibody_database](), true);
	start_database_session(default_path);

}

AntibodyDatabaseManager::AntibodyDatabaseManager(std::string const database_path) {
	start_database_session(database_path);
}

AntibodyDatabaseManager::~AntibodyDatabaseManager() {}
    
void
AntibodyDatabaseManager::start_database_session(std::string const database_path) {
	db_path_ = database_path;
	db_session_ = basic::database::get_db_session(db_path_);         
}

vector1< CDRNameEnum >
AntibodyDatabaseManager::load_cdr_design_data(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, std::map<core::Size,AAProbabilities>& prob_set, core::Size const cutoff, DesignInstructions & instructions) {
	if (! ab_info->clusters_setup()){
		utility_exit_with_message("Cluster information must be set in AntibodyInfo to load the design probabilities");
	}
	
	TR << "Loading CDR cluster statistics " << std::endl;
	vector1<CDRNameEnum> cdrs_with_no_data;
	
	for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		
		if (! instructions[cdr].design) {
			continue;
		}
		else if (instructions[cdr].design && instructions[cdr].conservative_design){
			continue;
		}


		std::pair<CDRClusterEnum, core::Real > cluster_pair = ab_info->get_CDR_cluster(cdr);
		CDRClusterEnum cluster = cluster_pair.first;
		if (cluster == NA){
			cdrs_with_no_data.push_back(cdr);
			TR << ab_info->get_CDR_Name(cdr) << " is of unknown cluster.  No design probability data added. Using conservative mutations instead." << std::endl;
			continue;
		}
		
		//Check to make sure cluster is present in antibody database
		std::string base_statement =
			"SELECT \n"
			"	position, \n"
			"	probability, \n"
			"	aa, \n"
			"	total_seq \n"
			"FROM\n"
			"	cdr_residue_probabilities \n"
			"WHERE \n"
			"	fullcluster = ?";
		
		cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
		select_statement.bind(1, ab_info->get_cluster_name(cluster));
		cppdb::result prob_result(basic::database::safely_read_from_database(select_statement));
		

		core::Size total_seq;
		while(prob_result.next()){
			
			if (prob_result.empty()) {
				TR << ab_info->get_CDR_Name(cdr) << " has no design probability data.  Using conservative mutations instead." << std::endl;
				break;
			}
			
			core::Size position;
			core::Real probability;
			std::string amino;
			prob_result >> position >> probability >> amino >> total_seq;
			
			
			if (  total_seq < cutoff){
				cdrs_with_no_data.push_back(cdr);
				TR << " Total data points for " << ab_info->get_cluster_name(cluster)<<" at "<< total_seq << " is lower than the set cutoff value of "<< cutoff << std::endl;
				TR << " Using conservative mutations instead." << std::endl;
				break;
			}
			
			
			core::Size rosetta_resnum = ab_info->get_CDR_start(cdr, pose) + position - 1;
			prob_set[rosetta_resnum][core::chemical::aa_from_oneletter_code(amino[0])] = probability;
		}
		TR << "Loaded "<< ab_info->get_cluster_name(cluster) << " with " << total_seq << " datapoints. " << std::endl;
	}
	return cdrs_with_no_data;
}

void
AntibodyDatabaseManager::check_for_graft_instruction_inconsistencies(AntibodyInfoCOP ab_info, GraftInstructions & instructions) {
	//Double check to make sure everything is kosher before grabbing Poses.
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if ((instructions[cdr].include_only_clusters.size() >= 1) && (instructions[cdr].leave_out_clusters.size() >= 1)){
			for(core::Size j=1; j<=instructions[cdr].include_only_clusters.size(); ++j){
				for (core::Size k=1; k<=instructions[cdr].leave_out_clusters.size(); ++k){
					CDRClusterEnum include_cluster = instructions[cdr].include_only_clusters[j];
					CDRClusterEnum leave_cluster = instructions[cdr].leave_out_clusters[k];
					
					//ab_info->get_cluster_name(include_cluster)
					if(include_cluster == leave_cluster){
						utility_exit_with_message("Cannot leave out and include cluster: "+ab_info->get_cluster_name(include_cluster));
					}
				}
			}
		}
		if ((instructions[cdr].include_only_pdb_ids.size() >= 1) && (instructions[cdr].leave_out_pdb_ids.size() >= 1)){
			for(core::Size j=1; j<=instructions[cdr].include_only_pdb_ids.size(); ++j){
				for (core::Size k=1; k<=instructions[cdr].leave_out_pdb_ids.size(); ++k){
					std::string include_pdbid = instructions[cdr].include_only_pdb_ids[j];
					std::string leave_pdbid = instructions[cdr].leave_out_pdb_ids[k];
					
					if (include_pdbid == leave_pdbid){
						utility_exit_with_message("Cannot leave out and include PDB id: "+include_pdbid);
					}
				}
			}
		}
	}
}

std::pair<CDRSet, CDRClusterMap>
AntibodyDatabaseManager::load_cdrs_for_grafting(AntibodyInfoCOP ab_info, GraftInstructions& instructions, PDBMap & pdbmap, core::Size overhang /* 3 */){
	
	
	//Check to make sure everything is kosher
	check_for_graft_instruction_inconsistencies(ab_info, instructions);
	
	
	protocols::features::ProteinSilentReport reporter = protocols::features::ProteinSilentReport();
	CDRSet cdr_set;
	CDRClusterMap cdr_cluster_map;
	
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		// Skip Fixed CDRs.  Yay for maps.
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (!instructions[cdr].graft){continue;}
		
		std::string base_statement = 
			"SELECT\n"
			"	PDB,\n "
			"	fullcluster,\n"
			"	original_chain\n"
			"FROM\n"
			"	cdr_data \n"
			"WHERE\n"
			"	CDR = ? AND\n"
			"	center = ? AND\n"
			"	length >= ? AND\n"
			"	length <= ? AND\n"
			//"	cluster != ? AND\n"
			"	datatag != ?";
		
		//Optional part of Query
		if (instructions[cdr].stay_native_cluster){
			
			//Remove incompatible settings.
			instructions[cdr].cluster_centers_only = false;
			instructions[cdr].leave_out_clusters.clear();
			instructions[cdr].include_only_clusters.clear();
			for (core::Size type = 1; type <=3 ; ++type){
				instructions[cdr].cluster_types[type]=true;
			}
			
			TR << ab_info->get_CDR_Name(cdr)<< 
					"Staying native cluster.  Ignoring CLUSTER + CENTER settings for this CDR.\n" <<
					"If PDBID INCLUDEONLY is set, please check that these PDBs have CDRs that match your cluster\n"<<std::endl;
			base_statement += " AND fullcluster=?\n";
				
		}
		
		//H3 Currently has NO types.
		if (cdr!=h3){
			base_statement += " AND type IN";
			
			//I don't know how to do this more pretty
			
			core::Size total_true=0;
			for (core::Size type = 1; type<=3; ++type){
				if (instructions[cdr].cluster_types[type]){
					total_true+=1;
				}
			}
			base_statement += get_string_for_IN(total_true);
		}
		
		//Optional part of Query
		//Using the IN keyword, could have just as easily done OR
		if (instructions[cdr].include_only_clusters.size() >= 1){
			base_statement += " AND fullcluster IN";
			base_statement += get_string_for_IN(instructions[cdr].include_only_clusters.size());
		}
		if (instructions[cdr].leave_out_clusters.size() >= 1){
			base_statement += " AND fullcluster NOT IN";
			base_statement += get_string_for_IN(instructions[cdr].leave_out_clusters.size());
		}
		if (instructions[cdr].include_only_pdb_ids.size() >= 1){
			base_statement += " AND PDB IN";
			base_statement += get_string_for_IN(instructions[cdr].include_only_pdb_ids.size());
		}
		if (instructions[cdr].leave_out_pdb_ids.size() >= 1){
			base_statement += " AND PDB NOT IN";
			base_statement += get_string_for_IN(instructions[cdr].leave_out_pdb_ids.size());
		}
		
		//Add stuff from Instructions
		base_statement = base_statement+";";
		TR.Debug<<"Final Statement: "<< base_statement << std::endl;
		cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
		select_statement.bind(1, ab_info->get_CDR_Name(cdr));
		select_statement.bind(2, (int)instructions[cdr].cluster_centers_only);
		select_statement.bind(3, instructions[cdr].min_length); //This SHOULD get compared to int.
		select_statement.bind(4, instructions[cdr].max_length);
		//select_statement.bind(5, -1);
		select_statement.bind(5, "loopKeyNotInPaper");
		
		//Optional Binds + Cluster types
		core::Size col = 5; //Keeps track of where we are.
		if (cdr != h3){
			for (core::Size t=1; t <= 3; ++t){
				if (instructions[cdr].cluster_types[t]){
					col+=1;
					select_statement.bind(col, t);
				}
			}
		}
		if (instructions[cdr].stay_native_cluster){
			col += 1;
			std::pair<CDRClusterEnum, core::Real> cluster = ab_info->get_CDR_cluster(cdr);
			if(cluster.first == NA){
				utility_exit_with_message(ab_info->get_CDR_Name(cdr)+" : Unable to identify cluster.  Modify StayNativeCluster in instructions or download new AntibodyDatabase + update Rosetta for new definitions.");
			select_statement.bind(col, ab_info->get_cluster_name(cluster.first));
			}
		}
		
		if (instructions[cdr].include_only_clusters.size() >= 1){
			for (core::Size j = 1; j <= instructions[cdr].include_only_clusters.size(); ++j){
				col += 1;
				select_statement.bind(col, ab_info->get_cluster_name(instructions[cdr].include_only_clusters[j]));
			}
		}
		if (instructions[cdr].leave_out_clusters.size() >= 1){
			for (core::Size j = 1; j <= instructions[cdr].leave_out_clusters.size(); ++j){
				col += 1;
				select_statement.bind(col, ab_info->get_cluster_name(instructions[cdr].leave_out_clusters[j]));
			}
		}
		if (instructions[cdr].include_only_pdb_ids.size() >= 1){
			for (core::Size j = 1; j <= instructions[cdr].include_only_pdb_ids.size(); ++j){
				col += 1;
				select_statement.bind(col, instructions[cdr].include_only_pdb_ids[j]);
			}
		}
		if (instructions[cdr].leave_out_pdb_ids.size() >=1 ){
			for (core::Size j = 1; j <= instructions[cdr].leave_out_pdb_ids.size(); ++j){
				col += 1;
				select_statement.bind(col, instructions[cdr].leave_out_pdb_ids[j]);
			}
		}
		
		
		//Get Result + Match to Rosetta structure database.  Create PoseOPs for the CDRSet.
		cppdb::result pdbid_result(basic::database::safely_read_from_database(select_statement));
		vector1<std::string> tags;
		vector1<CDRClusterEnum > clusters;
		CDRClusterEnumManagerOP manager = ab_info->get_cdr_cluster_enum_manager();
		while(pdbid_result.next()){
			std::string pdbid;
			std::string cluster;
			std::string original_chain;
			pdbid_result >> pdbid >> cluster >> original_chain;
			
			//Convert from my database to Rosetta input
			boost::to_lower(pdbid); //Original structures are lower case to distinguish chain easily in filename.
			std::string pdbid_in = pdbid+original_chain+"_"+ab_info->get_CDR_Name(cdr)+"_0001";
			
			if (manager->cdr_cluster_is_present(cluster)){
				
				CDRClusterEnum cluster_enum = ab_info->get_cluster_enum(cluster);
				clusters.push_back(cluster_enum);
				tags.push_back(pdbid_in);
				TR.Debug << "PDBId: " << pdbid<<std::endl;
			}
			else{
				TR<< cluster << " Not Present in enums! Skipping" << std::endl;
			}

		}
		
		TR.Debug << "Total Poses: "<< tags.size()<<std::endl;
		//Get a single struct_id
		for (core::Size k=1; k<=tags.size(); ++k){
			std::string statement = 
				"SELECT\n"
				"	struct_id\n"
				"FROM\n"
				"	structures \n"
				"WHERE\n"
				"	tag = ?;";
			
			cppdb::statement select_statement(basic::database::safely_prepare_statement(statement, db_session_));
			select_statement.bind(1, tags[k]);
			cppdb::result uuid_result(basic::database::safely_read_from_database(select_statement));
			//if (uuid_result.empty()){continue;}
			std::string tag = tags[k];
			//Should only be one result.  But, this should protect us from non-matching info between data and structures for whatever reason.
			while(uuid_result.next()){
				if (uuid_result.empty()){break;}
				TR.Debug << "Loading structure: " << tags[k] <<std::endl;
				StructureID struct_id;
				uuid_result >> struct_id;
				core::pose::PoseOP pose = new core::pose::Pose();
				reporter.load_pose(db_session_, struct_id, *pose);
				pose->conformation().detect_disulfides();
				
				if (pose->total_residue() != ab_info->get_cluster_length(clusters[k])+ overhang*2 ) {
					TR << "Leaving out bad structure:  "<< tags[k] << " num residues in pose do not match cluster cdr length with overhang. " <<std::endl;
					continue;
				}
				cdr_set[cdr].push_back(pose);
				cdr_cluster_map[cdr].push_back(clusters[k]);
				pdbmap[cdr].push_back(tag);
			}
		}
	}
	std::pair<CDRSet, CDRClusterMap> set_and_map = std::make_pair(cdr_set, cdr_cluster_map);
	return set_and_map;
}

 

} //design
} //antibody
} //protocols
