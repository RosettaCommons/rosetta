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
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>

//Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AA.hh>

//Database Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <protocols/features/ProteinSilentReport.hh>
#include <protocols/features/util.hh>
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
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/py/PyAssert.hh>
#include <boost/algorithm/string.hpp>
#include <fstream>


static thread_local basic::Tracer TR("antibody.design.AntibodyDatabaseManager");

namespace protocols {
namespace antibody {

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
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::design;

typedef std::map< CDRNameEnum, vector1<CDRPose> > CDRSet;


AntibodyDatabaseManager::AntibodyDatabaseManager(AntibodyInfoCOP ab_info, bool force_north_paper_db):
	utility::pointer::ReferenceCount(),
	ab_info_(ab_info)
{
	//protein_silent_report_ = new protocols::features::ProteinSilentReport();
	//ab_info_ = ab_info;

	string default_path = basic::database::full_name(basic::options::option[basic::options::OptionKeys::antibody::design::antibody_database](), true);
	string north_paper_path = basic::database::full_name(basic::options::option[basic::options::OptionKeys::antibody::design::paper_ab_db_path](), true);
	bool use_north_paper_ab_db = basic::options::option[basic::options::OptionKeys::antibody::design::paper_ab_db]();

	if (use_north_paper_ab_db || force_north_paper_db){
		start_database_session(north_paper_path);
	}
	else if (utility::file::file_exists(default_path)){
		start_database_session(default_path);
	}
	else {
		utility_exit_with_message("Could not locate antibody database.  Please check paths, use the -paper_ab_db option, or pass the class constructor option to force the use of the paper database.");
	}


}

AntibodyDatabaseManager::AntibodyDatabaseManager(AntibodyInfoCOP ab_info, std::string const database_path):
	utility::pointer::ReferenceCount(),
	ab_info_(ab_info)
{
	//ab_info_ = ab_info;
	start_database_session(database_path);
}

AntibodyDatabaseManager::~AntibodyDatabaseManager() {}

void
AntibodyDatabaseManager::start_database_session(std::string const database_path) {
	db_path_ = database_path;
	db_session_ = basic::database::get_db_session(db_path_);
}

vector1< CDRNameEnum >
AntibodyDatabaseManager::load_cdr_design_data(AntibodyCDRSeqDesignOptions const & instructions, core::pose::Pose const & pose, std::map<core::Size,AAProbabilities>& prob_set, core::Size const cutoff) {
	TR << "Loading CDR cluster statistics " << std::endl;
	vector1<CDRNameEnum> cdrs_with_no_data;

	for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		CDRSeqDesignOptionsCOP options = instructions[cdr];

		if (! options->design()) {
			continue;
		}
		else if (options->design() && options->design_strategy() != seq_design_profiles){
			continue;
		}


		// Load data from either the stored DataCache or AntibodyInfo.
		CDRClusterEnum cluster;

		if (pose.data().has(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO)){
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
			cluster = cluster_cache.get_cluster(cdr)->cluster();
		}
		else {
			cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
		}

		if (cluster == NA ){
			cdrs_with_no_data.push_back(cdr);
			TR << ab_info_->get_CDR_name(cdr) << " is of unknown cluster.  No design probability data added. Using conservative mutations instead." << std::endl;
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
		select_statement.bind(1, ab_info_->get_cluster_name(cluster));
		cppdb::result prob_result(basic::database::safely_read_from_database(select_statement));


		core::Size total_seq;
		while(prob_result.next()){

			if (prob_result.empty()) {
				TR << ab_info_->get_CDR_name(cdr) << " has no design probability data.  Using conservative mutations instead." << std::endl;
				break;
			}

			core::Size position;
			core::Real probability;
			std::string amino;
			prob_result >> position >> probability >> amino >> total_seq;


			if (  total_seq < cutoff){
				cdrs_with_no_data.push_back(cdr);
				TR << " Total data points for " << ab_info_->get_cluster_name(cluster)<<" at "<< total_seq << " is lower than the set cutoff value of "<< cutoff << std::endl;
				TR << " Using conservative mutations instead." << std::endl;
				break;
			}


			core::Size rosetta_resnum = ab_info_->get_CDR_start(cdr, pose) + position - 1;
			prob_set[rosetta_resnum][core::chemical::aa_from_oneletter_code(amino[0])] = probability;
		}
		TR << "Loaded "<< ab_info_->get_cluster_name(cluster) << " with " << total_seq << " datapoints. " << std::endl;
	}
	return cdrs_with_no_data;
}

void
AntibodyDatabaseManager::check_for_graft_instruction_inconsistencies(AntibodyCDRSetOptions const & instructions) {

	//Double check to make sure everything is kosher before grabbing Poses.
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRSetOptionsCOP options = instructions[cdr];
		if (has_vec_inconsistency(options->include_only_clusters(), options->exclude_clusters())){
			utility_exit_with_message("Inconsistency in CLUSTERs option. Cannot leave out and include option.");
		}

		if (has_vec_inconsistency(options->include_only_pdbs(), options->exclude_pdbs())){
			utility_exit_with_message("Inconsistency in PDBID option. Cannot leave out and include option.");
		}

		if (has_vec_inconsistency(options->include_only_species(), options->exclude_species())){
			utility_exit_with_message("Inconsistency in SPECIES option. Cannot leave out and include option.");
		}

		if (has_vec_inconsistency(options->include_only_germlines(), options->exclude_germlines() )){
			utility_exit_with_message("Inconsistency in GERMLINE option. Cannot leave out and include option.");
		}
	}
}

//CDRSet
//AntibodyDatabaseManager::load_cdr_poses(const CDRSetOptions& options, const core::pose::Pose& pose, const bool use_light_chain_type = true, core::Size overhang = 3) {

//}

CDRSet
AntibodyDatabaseManager::load_cdr_poses(AntibodyCDRSetOptions const & instructions, core::pose::Pose const & pose, bool const use_light_chain_type, core::Size overhang /* 3 */){


	//Check to make sure everything is kosher
	check_for_graft_instruction_inconsistencies(instructions);

	TR << "loading cdrs from database " << db_path_  << std::endl;
	protocols::features::ProteinSilentReport reporter = protocols::features::ProteinSilentReport();
	CDRSet cdr_set;

	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		// Skip Fixed CDRs.  Yay for maps.
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRSetOptionsOP options = instructions[cdr]->clone();

		if (!options->load()){continue;}

		std::string base_statement =
			"SELECT\n"
			"	PDB,\n "
			"	fullcluster,\n"
			"	original_chain,\n"
			"	dis,\n"
			"	resolution\n"
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
		if (options->include_only_current_cluster()){

			//Remove incompatible settings.
			utility::vector1<CDRClusterEnum> empty_vector1;
			utility::vector1<CDRClusterEnum> empty_vector2;

			options->include_only_center_clusters(false);
			options->exclude_clusters(empty_vector1);
			options->include_only_clusters(empty_vector2);
			for (core::Size type = 1; type <=3 ; ++type){
				options->length_type(type, true);
			}

			TR << ab_info_->get_CDR_name(cdr)<<
					" Staying native cluster.  Ignoring CLUSTER + CENTER settings for this CDR.\n" <<
					"If PDBID INCLUDEONLY is set, please check that these PDBs have CDRs that match your cluster\n"<<std::endl;
			base_statement += " AND fullcluster=?\n";

		}

		//H3 Currently has NO types.
		if (cdr!=h3){
			base_statement += " AND length_type IN";

			//I don't know how to do this more pretty

			core::Size total_true=0;
			for (core::Size type = 1; type<=3; ++type){
				if (options->length_type()[type]){
					total_true+=1;
				}
			}
			base_statement += get_question_mark_string(total_true);
		}

		//Optional part of Query - Not sure if better way than a bunch of if statements
		//Using the IN keyword, could have just as easily done OR

		//Cluster
		if (options->include_only_clusters().size() >= 1){
			base_statement += " AND fullcluster IN";
			base_statement += get_question_mark_string(options->include_only_clusters().size());
		}
		if (options->exclude_clusters().size() >= 1){
			base_statement += " AND fullcluster NOT IN";
			base_statement += get_question_mark_string(options->exclude_clusters().size());
		}

		//PDBID
		if (options->include_only_pdbs().size() >= 1){
			base_statement += " AND PDB IN";
			base_statement += get_question_mark_string(options->include_only_pdbs().size());
		}
		if (options->exclude_pdbs().size() >= 1){
			base_statement += " AND PDB NOT IN";
			base_statement += get_question_mark_string(options->exclude_pdbs().size());
		}

		//Species
		if (options->include_only_species().size() >= 1){
			base_statement += " AND GSpecies IN"; //IMGT Germline determined CDR species - two letter code in docs
			base_statement += get_question_mark_string(options->include_only_species().size());
		}
		if (options->exclude_species().size() >= 1) {
			base_statement += " AND GSpecies NOT IN ";
			base_statement += get_question_mark_string(options->exclude_species().size());
		}

		//GermLine
		if (options->include_only_germlines().size() >= 1) {
			base_statement += " AND IG IN";
			base_statement += get_question_mark_string(options->include_only_germlines().size());
		}
		if (options->exclude_germlines().size() >= 1){
			base_statement += " AND IG IN";
			base_statement += get_question_mark_string(options->exclude_germlines().size());
		}
		//Lambda vs Kappa
		if (ab_info_->get_light_chain_type_enum() != unknown && use_light_chain_type){
			base_statement += " AND (gene = 'heavy' OR gene = '"+ab_info_->get_light_chain_type()+"')";
		}

		//Add stuff from Instructions
		base_statement = base_statement+";";
		if (TR.Debug.visible()){
			TR.Debug << "Final Statement:\n"<< base_statement << std::endl;
		}
		cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
		select_statement.bind(1, ab_info_->get_CDR_name(cdr));
		select_statement.bind(2, (int)options->include_only_center_clusters());
		select_statement.bind(3, options->min_length()); //This SHOULD get compared to int.
		select_statement.bind(4, options->max_length());
		//select_statement.bind(5, -1);
		select_statement.bind(5, "loopKeyNotInPaper");

		//Optional Binds + Cluster types
		core::Size col = 5; //Keeps track of where we are.
		if (cdr != h3){
			for (core::Size t=1; t <= 3; ++t){
				if (options->length_type()[t]){
					col+=1;
					select_statement.bind(col, t);
				}
			}
		}
		if (options->include_only_current_cluster()){
			col += 1;

			// Load data from either the stored DataCache or AntibodyInfo.
			CDRClusterEnum current_cluster;

			if (pose.data().has(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO)){
				BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
				current_cluster = cluster_cache.get_cluster(cdr)->cluster();
			}
			else {
				current_cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
			}
			if(current_cluster  == NA){
				utility_exit_with_message(ab_info_->get_CDR_name(cdr)+" : Unable to identify cluster.  Modify StayNativeCluster in instructions or download new AntibodyDatabase + update Rosetta for new definitions.");
				select_statement.bind(col, ab_info_->get_cluster_name(current_cluster));
			}
		}

		//Bind optional string vector constraints

		bind_vec_constraint(get_cluster_string_vec(options->include_only_clusters()), select_statement, col);
		bind_vec_constraint(get_cluster_string_vec(options->exclude_clusters()), select_statement, col);

		bind_vec_constraint(options->include_only_pdbs(), select_statement, col);
		bind_vec_constraint(options->exclude_pdbs(), select_statement, col);

		bind_vec_constraint(options->include_only_species(), select_statement, col);
		bind_vec_constraint(options->exclude_species(), select_statement, col);

		bind_vec_constraint(options->include_only_germlines(), select_statement, col);
		bind_vec_constraint(options->exclude_germlines(), select_statement, col);

		//Get Result + Match to Rosetta structure database.  Create PoseOPs for the CDRSet.
		cppdb::result pdbid_result(basic::database::safely_read_from_database(select_statement));

		vector1<std::string> tags;
		vector1<CDRClusterEnum > clusters;
		vector1< core::Real > distances;
		vector1< core::Real > resolutions;

		CDRClusterEnumManagerCOP manager = ab_info_->get_cdr_cluster_enum_manager();
		while(pdbid_result.next()){
			std::string pdbid;
			std::string cluster;
			std::string original_chain;
			core::Real dis;
			core::Real resolution;

			pdbid_result >> pdbid >> cluster >> original_chain >> dis >> resolution;

			//TR << pdbid<<" "<<cluster << original_chain << std::endl;
			//Convert from my database to Rosetta input
			boost::to_lower(pdbid); //Original structures are lower case to distinguish chain easily in filename.
			std::string pdbid_in = pdbid+original_chain+"_"+ab_info_->get_CDR_name(cdr)+"_0001";

			if (manager->cdr_cluster_is_present(cluster)){

				CDRClusterEnum cluster_enum = ab_info_->get_cluster_enum(cluster);
				clusters.push_back(cluster_enum);

				//TR.Debug << "PDBId: " << pdbid<<std::endl;
			}
			else if (cdr== h3){
				//Increased H3 sampling
				clusters.push_back(NA);
				//TR.Debug << "PDBId: " << pdbid<<std::endl;

			}
			else {
				TR<< cluster << " Not Present in enums! Skipping" << std::endl;
				continue;
			}
			tags.push_back(pdbid_in);
			distances.push_back(dis);
			resolutions.push_back(resolution);

		}

		TR << "Total "<<ab_info_->get_CDR_name(cdr) << " poses: " << tags.size() << std::endl;

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
				//TR.Debug << "Loading structure: " << tags[k] <<std::endl;
				StructureID struct_id;
				uuid_result >> struct_id;
				core::pose::PoseOP pose( new core::pose::Pose() );
				reporter.load_pose(db_session_, struct_id, *pose);
				pose->conformation().detect_disulfides();

				if (pose->total_residue() != ab_info_->get_cluster_length(clusters[k])+ overhang*2 ) {
					TR << "Leaving out bad structure:  "<< tags[k] << " num residues in pose do not match cluster cdr length with overhang. " <<std::endl;
					continue;
				}
				CDRPose cdr_pose = CDRPose(pose, clusters[k], tag, distances[k], resolutions[k]);
				cdr_set[cdr].push_back(cdr_pose);

			}
		}
	}
	TR << "cdrs loaded from database."<<std::endl;
	return cdr_set;

	//std::pair<CDRSet, CDRClusterMap>
}

template < typename T >
bool
AntibodyDatabaseManager::has_vec_inconsistency(vector1<T> const &  include, vector1<T> const & exclude) const  {

	if ((include.size() >= 1) && (exclude.size() >= 1)){
		for(core::Size j=1; j<=include.size(); ++j){
			for (core::Size k=1; k<=exclude.size(); ++k){

				if (include[j]== exclude[k]){
					return true;
				}
			}
		}
	}

	return false;
}

void
AntibodyDatabaseManager::bind_vec_constraint(vector1<std::string> const & vec, cppdb::statement& select_statement, core::Size& col) const {

	//Could not get explicit template function specialization to
	if (vec.size() >= 1) {
		for (core::Size j = 1; j <= vec.size(); ++j){
			col +=1;
			select_statement.bind(col, vec[j]);
		}
	}
}

vector1<std::string>
AntibodyDatabaseManager::get_cluster_string_vec(const vector1<CDRClusterEnum>& clusters){
	vector1<std::string> clusters_str;
	for (core::Size i = 1; i <= clusters.size(); ++i){
		clusters_str.push_back(ab_info_->get_cluster_name(clusters[i]));
	}
	return clusters_str;
}

} //antibody
} //protocols
