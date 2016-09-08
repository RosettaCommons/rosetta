// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyDatabaseManager.cc
/// @brief Handles all loading of CDR, Framework, and cluster/dmap info from external database file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Antibody Headers
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyDesignMover.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>

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


static THREAD_LOCAL basic::Tracer TR("antibody.database.AntibodyDatabaseManager");

namespace protocols {
namespace antibody {

using std::string;
using std::endl;
using core::pose::Pose;
using utility::file::FileName;
using utility::vector1;
using utility::sql_database::DatabaseSessionManager;
using protocols::features::ProteinSilentReport;
using namespace core;
using namespace protocols::features;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::design;

typedef std::map< CDRNameEnum, vector1<CDRDBPose> > CDRDBPoseSet;


AntibodyDatabaseManager::AntibodyDatabaseManager(AntibodyInfoCOP ab_info, bool force_north_paper_db):
	utility::pointer::ReferenceCount(),
	ab_info_(ab_info)
{
	//protein_silent_report_ = new protocols::features::ProteinSilentReport();
	//ab_info_ = ab_info;

	string default_path = basic::database::full_name(basic::options::option[basic::options::OptionKeys::antibody::design::antibody_database](), true);
	string north_paper_path = basic::database::full_name(basic::options::option[basic::options::OptionKeys::antibody::design::paper_ab_db_path](), true);
	bool use_north_paper_ab_db = basic::options::option[basic::options::OptionKeys::antibody::design::paper_ab_db]();
	use_outliers_ = basic::options::option[basic::options::OptionKeys::antibody::design::use_outliers]();
	use_h3_graft_outliers_ = basic::options::option[basic::options::OptionKeys::antibody::design::use_H3_graft_outliers]();
	use_only_H3_kinked_ = basic::options::option[ basic::options::OptionKeys::antibody::design::use_only_H3_kinked]();
	use_light_chain_type_ = option [OptionKeys::antibody::design::use_light_chain_type]();

	if ( use_north_paper_ab_db || force_north_paper_db ) {
		start_database_session(north_paper_path);
	} else if ( utility::file::file_exists(default_path) ) {
		start_database_session(default_path);

	} else if ( utility::file::file_exists(north_paper_path) ) {
		start_database_session(north_paper_path);
	} else {
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

AntibodyDatabaseManager::~AntibodyDatabaseManager(){}

void
AntibodyDatabaseManager::start_database_session(std::string const database_path) {
	db_path_ = database_path;
	db_session_ = basic::database::get_db_session(db_path_);
	TR << "Reading from: " << db_path_ << std::endl;
}

void
AntibodyDatabaseManager::set_outlier_use(bool use_outliers) {
	use_outliers_ = use_outliers;
}

void
AntibodyDatabaseManager::use_light_chain_type(bool use_light_chain_type){
	use_light_chain_type_ = use_light_chain_type;
}

utility::vector1<core::Size>
AntibodyDatabaseManager::get_cluster_totals() const {

	utility::vector1<core::Size> totals = get_cluster_totals(use_outliers_);

	if ( use_h3_graft_outliers_ && !use_outliers_ ) {
		utility::vector1<core::Size> totals_w_outliers = get_cluster_totals(true);
		for ( core::Size i = 1; i <= core::Size(CDRClusterEnum_stop); ++i ) {
			CDRClusterEnum cluster = static_cast<CDRClusterEnum>( i );
			if ( ab_info_->get_cdr_enum_for_cluster(cluster) == h3 ) {
				totals[ cluster ] = totals_w_outliers[ cluster ];
			}
		}
	}
	return totals;
}

std::string
AntibodyDatabaseManager::get_ands(std::string name, core::Size n) const {
	if ( n == 1 ) {
		std::string out = name+"=?";
		return out;
	} else {
		std::string out = "("+name+"=?";
		for ( core::Size i = 2; i <= n; ++i ) {
			out = out +" OR "+name+"=?";
		}
		out = out+")";
		return out;
	}

}

utility::vector1<core::Size>
AntibodyDatabaseManager::get_cluster_totals(bool use_outliers) const{
	utility::vector1<core::Size> totals(CDRClusterEnum_total, 0);

	CDRClusterEnumManagerCOP manager = ab_info_->get_cdr_cluster_enum_manager();

	std::string seq_table = use_outliers ? "cdr_res_prob_outliers_true" : "cdr_res_prob_outliers_false_liberal";
	//Check to make sure cluster is present in antibody database
	std::string base_statement =
		"SELECT \n"
		"\tfullcluster, \n"
		"\ttotal_seq \n"
		"FROM\n"
		+ seq_table+"\n"
		"order by fullcluster";

	//TR << "Reading from: " << db_path_ << std::endl;

	cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
	cppdb::result total_result(basic::database::safely_read_from_database(select_statement));


	//core::Size total_seq;
	while ( total_result.next() ) {

		if ( total_result.empty() ) {
			break;
		}

		core::Size total_seq;
		std::string fullcluster;
		total_result >> fullcluster >> total_seq;

		if ( manager->cdr_cluster_is_present(fullcluster) ) {

			CDRClusterEnum cluster = ab_info_->get_cluster_enum(fullcluster);
			totals[ cluster ] = total_seq;

		} else {
			utility_exit_with_message("Cluster "+fullcluster+" not present in Enums!!");
		}
	}

	///H3 Outlier change here.
	return totals;

}

CDRDBPoseSet
AntibodyDatabaseManager::load_cdr_poses(
	AntibodyCDRSetOptions const & instructions,
	core::pose::Pose const & pose,
	core::Size overhang /* 3 */){


	//Check to make sure everything is kosher
	check_for_graft_instruction_inconsistencies(instructions);
	protocols::features::ProteinSilentReport reporter = protocols::features::ProteinSilentReport();

	utility::vector1<core::Size> totals = this->get_cluster_totals();

	CDRDBPoseSet cdr_set;

	for ( core::Size i=1; i<=CDRNameEnum_total; ++i ) {

		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRSetOptionsOP options = instructions[cdr]->clone();

		if ( !options->load() ) {
			continue;
		}

		std::string base_statement =
			"SELECT\n"
			"\tPDB,\n "
			"\tfullcluster,\n"
			"\toriginal_chain,\n"
			"\tdis,\n"
			"\tDistDegree,\n"
			"\tresolution,\n"
			"\trama\n"
			"FROM\n"
			"\tcdr_data \n"
			"WHERE\n"
			"\tCDR = ? AND\n"
			"\tlength >= ? AND\n"
			"\tlength <= ? AND\n"
			//" cluster != ? AND\n"
			"\tdatatag != ?";


		//Optional part of Query
		bool use_outliers;
		if ( cdr == h3 ) {
			use_outliers = use_h3_graft_outliers_;
		} else {
			use_outliers = use_outliers_;
		}

		if ( ! use_outliers && ! options->include_only_center_clusters() ) {
			base_statement = base_statement +"\n"+
				"          AND (bb_rmsd_cdr_align < 1.5 OR DistDegree < 40)\n"
				"          AND DistDegree != -1\n"
				"          AND bb_rmsd_cdr_align != -1\n";
		}


		if ( options->include_only_current_cluster() ) {

			TR << ab_info_->get_CDR_name(cdr)<<
				" Staying native cluster.  Ignoring CENTER, CLUSTER, LENGTH, LENGTH_TYPE, and CLUSTER_CUTOFF settings for this CDR.\n" <<
				"If PDBID INCLUDEONLY is set, please check that these PDBs have CDRs that match your cluster\n"<<std::endl;

			//Remove incompatible settings.
			options->include_only_center_clusters(false);
			options->max_length(100);
			options->min_length(1);

			options->exclude_clusters_clear();
			options->include_only_clusters_clear();
			options->cluster_sampling_cutoff( 0 );

			for ( core::Size type = 1; type <=3 ; ++type ) {
				options->length_type(type, true);
			}


			// Load data from either the stored DataCache or AntibodyInfo.
			CDRClusterEnum current_cluster;

			if ( pose.data().has(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO) ) {
				BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
				current_cluster = cluster_cache.get_cluster(cdr)->cluster();
			} else {
				current_cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
			}
			if ( current_cluster  == NA ) {
				utility_exit_with_message(ab_info_->get_CDR_name(cdr)+" : Unable to identify cluster.  Modify StayNativeCluster in instructions or download new AntibodyDatabase + update Rosetta for new definitions.");

			}
			TR <<"Staying native cluster: " << ab_info_->get_cluster_name(current_cluster) << std::endl;
			options->include_only_clusters_add(current_cluster);

		}

		//H3 Currently has NO types.
		if ( cdr!=h3 ) {
			base_statement += " AND length_type IN";

			//I don't know how to do this more pretty

			core::Size total_true=0;
			for ( core::Size type = 1; type<=3; ++type ) {
				if ( options->length_type()[type] ) {
					total_true+=1;
				}
			}
			base_statement += get_question_mark_string(total_true);
		}

		//Optional part of Query - Not sure if better way than a bunch of if statements
		//Using the IN keyword, could have just as easily done OR

		//Cluster
		if ( options->include_only_clusters().size() >= 1 ) {
			base_statement += " AND fullcluster IN";
			base_statement += get_question_mark_string(options->include_only_clusters().size());
		}
		if ( options->exclude_clusters().size() >= 1 ) {
			base_statement += " AND fullcluster NOT IN";
			base_statement += get_question_mark_string(options->exclude_clusters().size());
		}

		//PDBID
		if ( options->include_only_pdbs().size() >= 1 ) {
			base_statement += " AND PDB IN";
			base_statement += get_question_mark_string(options->include_only_pdbs().size());
			//base_statement += " AND ";
			//base_statement += this->get_ands("PDB", options->include_only_pdbs().size());
		}
		if ( options->exclude_pdbs().size() >= 1 ) {
			base_statement += " AND PDB NOT IN";
			base_statement += get_question_mark_string(options->exclude_pdbs().size());
		}

		//Species
		if ( options->include_only_species().size() >= 1 ) {
			base_statement += " AND GSpecies IN"; //IMGT Germline determined CDR species - two letter code in docs
			base_statement += get_question_mark_string(options->include_only_species().size());
		}
		if ( options->exclude_species().size() >= 1 ) {
			base_statement += " AND GSpecies NOT IN ";
			base_statement += get_question_mark_string(options->exclude_species().size());
		}

		//GermLine
		if ( options->include_only_germlines().size() >= 1 ) {
			base_statement += " AND IG IN";
			base_statement += get_question_mark_string(options->include_only_germlines().size());
		}
		if ( options->exclude_germlines().size() >= 1 ) {
			base_statement += " AND IG IN";
			base_statement += get_question_mark_string(options->exclude_germlines().size());
		}
		//Lambda vs Kappa
		if ( ab_info_->get_light_chain_type_enum() != unknown && use_light_chain_type_ && (ab_info_->get_light_chain_type_enum() != lambda) ) {
			base_statement += " AND (gene = 'heavy' OR gene = '"+ab_info_->get_light_chain_type()+"')";
		} else if ( ab_info_->get_light_chain_type_enum() == lambda && use_light_chain_type_ ) {
			base_statement += " AND (gene = 'heavy' OR gene = 'lambda6' OR gene = 'lambda')";
		}
		if ( options->include_only_center_clusters() ) {
			base_statement = base_statement +
				"         AND center = 1";
		}

		//Add stuff from Instructions
		base_statement = base_statement+";";

		//TR<< "Final Statement:\n"<< base_statement << std::endl;


		cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
		select_statement.bind(1, ab_info_->get_CDR_name(cdr));

		select_statement.bind(2, options->min_length()); //This SHOULD get compared to int.
		select_statement.bind(3, options->max_length());
		//select_statement.bind(5, -1);
		select_statement.bind(4, "loopKeyNotInPaper");

		//Optional Binds + Cluster types
		core::Size col = 4; //Keeps track of where we are.
		if ( cdr != h3 ) {
			for ( core::Size t=1; t <= 3; ++t ) {
				if ( options->length_type()[t] ) {
					col+=1;
					select_statement.bind(col, t);
				}
			}
		}

		//Bind optional string vector constraints

		//This needs to happen or binding will not occur - strange thing with binding - the variable cannot be returned
		// and put into sql in C++, needs to be stored somewhere first.  Thank you Rocco.

		utility::vector1<std::string> include_clusters = get_cluster_string_vec(options->include_only_clusters());
		utility::vector1<std::string> exclude_clusters = get_cluster_string_vec(options->exclude_clusters());

		//TR << "include clusters" << utility::to_string(include_clusters) << std::endl;
		//TR << "exclude clusters" << utility::to_string(exclude_clusters) << std::endl;

		bind_vec_constraint(include_clusters, select_statement, col);
		bind_vec_constraint(exclude_clusters, select_statement, col);

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
		vector1< core::Real > normalized_distances;
		vector1< core::Real > resolutions;

		core::Size non_kinked_removed =0;
		vector1< core::Size > clusters_removed(CDRClusterEnum_stop, 0);


		if ( options->cluster_sampling_cutoff() > 0 ) {
			TR << "Removing all CDRs with cluster size < " << options->cluster_sampling_cutoff() << std::endl;
		}

		CDRClusterEnumManagerCOP manager = ab_info_->get_cdr_cluster_enum_manager();

		while ( pdbid_result.next() ) {
			std::string pdbid;
			std::string cluster;
			std::string original_chain;
			core::Real dis;
			core::Real normalized_dis;
			core::Real resolution;
			std::string rama;

			pdbid_result >> pdbid >> cluster >> original_chain >> dis >> normalized_dis >> resolution >> rama;

			// Skip if not kinked and we only want kinked here!
			if ( cdr == h3 && use_only_H3_kinked_ && ! is_H3_rama_kinked(rama) ) {
				non_kinked_removed+=1;
				continue;
			}
			//TR << pdbid<<" "<<cluster << original_chain << std::endl;
			//Convert from my database to Rosetta input
			boost::to_lower(pdbid); //Original structures are lower case to distinguish chain easily in filename.
			std::string pdbid_in = pdbid+original_chain+"_"+ab_info_->get_CDR_name(cdr)+"_0001";

			if ( manager->cdr_cluster_is_present(cluster) ) {
				CDRClusterEnum cluster_enum = ab_info_->get_cluster_enum(cluster);
				if ( options->cluster_sampling_cutoff() != 0 ) {
					if ( totals[ cluster_enum ] < options->cluster_sampling_cutoff() ) {
						TR << "Removed "<< cluster << " at " << totals[ cluster_enum ] << std::endl;
						clusters_removed[ cluster_enum ] = clusters_removed[ cluster_enum ] + 1;
						continue;
					}
				}

				clusters.push_back(cluster_enum);

				//TR.Debug << "PDBId: " << pdbid<<std::endl;
			} else if ( cdr== h3 ) {
				//Increased H3 sampling
				clusters.push_back(NA);
				//TR.Debug << "PDBId: " << pdbid<<std::endl;

			} else {
				TR<< cluster << " Not Present in enums! Skipping" << std::endl;
				continue;
			}
			tags.push_back(pdbid_in);
			distances.push_back(dis);
			resolutions.push_back(resolution);
			normalized_distances.push_back(normalized_dis);

		}

		// Print totals and counts of removed structures and clusters.
		TR << "Total "<<ab_info_->get_CDR_name(cdr) << " poses: " << tags.size() << std::endl;
		if ( non_kinked_removed > 0 ) {
			TR << "Removed "<< non_kinked_removed << " non-kinked H3 structures" << std::endl;
		}
		if ( options->cluster_sampling_cutoff() != 0 ) {
			core::Size total_removed_clusters = 0;
			core::Size total_removed_sequences = 0;

			for ( core::Size i = 1; i <= clusters_removed.size(); ++i ) {
				if ( clusters_removed[ i ] > 0 ) {
					total_removed_clusters+=1;
					total_removed_sequences = total_removed_sequences + clusters_removed[ i ];
				}
			}
			if ( total_removed_clusters > 0 ) {
				TR << "Removed " << total_removed_clusters << " small clusters" << std::endl;
			}
			if ( total_removed_sequences > 0 ) {
				TR << "Removed " << total_removed_sequences << " small cluster structures from CDRSet. " << std::endl;
			}
		}

		//Get a single struct_id
		for ( core::Size k=1; k<=tags.size(); ++k ) {
			std::string statement =
				"SELECT\n"
				"\tstruct_id\n"
				"FROM\n"
				"\tstructures \n"
				"WHERE\n"
				"\ttag = ?;";

			cppdb::statement select_statement(basic::database::safely_prepare_statement( statement, db_session_ ));
			select_statement.bind(1, tags[ k ]);
			cppdb::result uuid_result(basic::database::safely_read_from_database( select_statement ));
			//if (uuid_result.empty()){continue;}
			std::string tag = tags[ k ];
			//Should only be one result.  But, this should protect us from non-matching info between data and structures for whatever reason.
			while ( uuid_result.next() ) {
				if ( uuid_result.empty() ) { break; }
				//TR.Debug << "Loading structure: " << tags[k] <<std::endl;
				StructureID struct_id;
				uuid_result >> struct_id;
				core::pose::PoseOP pose( new core::pose::Pose() );
				reporter.load_pose( db_session_, struct_id, *pose );
				pose->conformation().detect_disulfides();

				if ( pose->size() != ab_info_->get_cluster_length(clusters[ k ])+ overhang*2 ) {
					TR << "Leaving out bad structure:  "<< tags[k] << " num residues in pose do not match cluster cdr length with overhang. " <<std::endl;
					continue;
				}
				CDRDBPose cdr_pose = CDRDBPose(pose, clusters[ k ], tag, distances[ k ], normalized_distances[ k ], resolutions[ k ]);
				cdr_set[cdr].push_back(cdr_pose);

			}
		}
	}
	TR << "cdrs loaded from database."<<std::endl;
	return cdr_set;

	//std::pair<CDRDBPoseSet, CDRClusterMap>
}

vector1< bool >
AntibodyDatabaseManager::load_cdr_design_data(
	AntibodyCDRSeqDesignOptions const & instructions,
	core::pose::Pose const & pose,
	std::map<core::Size,AAProbabilities>& prob_set,
	core::Size const cutoff)
{
	TR << "Loading CDR cluster statistics " << std::endl;
	vector1<bool>cdrs_to_load(6, true);

	for ( core::Size i = 1; i<=CDRNameEnum_total; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		CDRSeqDesignOptionsCOP options = instructions[cdr];

		if ( ! options->design() ) {
			cdrs_to_load[ cdr ] = false;
		} else if ( options->design() && options->design_strategy() != seq_design_profiles ) {
			cdrs_to_load[ cdr ] = false;
		}
	}
	return this->load_cdr_design_data_for_cdrs(cdrs_to_load, pose, prob_set, cutoff);
}

vector1< bool >
AntibodyDatabaseManager::load_cdr_design_data_for_cdrs(
	utility::vector1<bool> const & c,
	const core::pose::Pose& pose,
	std::map<core::Size,AAProbabilities>& prob_set,
	const core::Size cutoff)
{

	utility::vector1< bool > cdrs = c;
	if ( cdrs.size() == 6 ) {
		cdrs.push_back(false);
		cdrs.push_back(false);
	}


	assert(cdrs.size() == 8);
	vector1<bool> cdrs_with_no_data(8, false);

	for ( core::Size i = 1; i <= CDRNameEnum_total; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if ( ! cdrs[ cdr ] || ! ab_info_->has_CDR(cdr) ) continue;

		// Load data from either the stored DataCache or AntibodyInfo.
		CDRClusterEnum cluster;

		if ( pose.data().has(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO) ) {
			//TR <<"Setting up profile data for cluster from pose datacache "<< ab_info_->get_CDR_name( cdr ) << std::endl;
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
			cluster = cluster_cache.get_cluster(cdr)->cluster();
		} else {
			cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
		}

		if ( cluster == NA ) {
			cdrs_with_no_data[ cdr] = true;
			TR << ab_info_->get_CDR_name(cdr) << " is of unknown cluster.  No design probability data added."<< std::endl;
			continue;
		}

		std::string seq_table = use_outliers_ ? "cdr_res_prob_outliers_true" : "cdr_res_prob_outliers_false_liberal";
		//Check to make sure cluster is present in antibody database
		std::string base_statement =
			"SELECT \n"
			"\tposition, \n"
			"\tprobability, \n"
			"\taa, \n"
			"\ttotal_seq \n"
			"FROM\n"
			+ seq_table+"\n"
			"WHERE \n"
			"\tfullcluster = ?";

		//TR << "Reading from: " << db_path_ << std::endl;

		cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
		select_statement.bind(1, ab_info_->get_cluster_name(cluster));
		cppdb::result prob_result(basic::database::safely_read_from_database(select_statement));


		core::Size total_seq;
		while ( prob_result.next() ) {

			if ( prob_result.empty() ) {
				TR << ab_info_->get_CDR_name(cdr) << " has no design probability data." << std::endl;
				break;
			}

			core::Size position;
			core::Real probability;
			std::string amino;
			prob_result >> position >> probability >> amino >> total_seq;


			if (  total_seq < cutoff ) {
				cdrs_with_no_data[ cdr ] = true;
				TR << "Skipped " << ab_info_->get_cluster_name(cluster) <<" with "<< total_seq << " data point  "<< std::endl;
				break;
			}


			core::Size rosetta_resnum = ab_info_->get_CDR_start(cdr, pose, North) + position - 1;
			prob_set[rosetta_resnum][core::chemical::aa_from_oneletter_code(amino[0])] = probability;
		}
		if (  total_seq >= cutoff ) {
			TR << "Loaded "<< ab_info_->get_cluster_name(cluster) << " with " << total_seq << " datapoints. " << std::endl;
		}
	}

	return cdrs_with_no_data;
}

CDRDBSequenceSet
AntibodyDatabaseManager::load_cdr_sequences(
	utility::vector1<bool> const & c,
	const core::pose::Pose& pose,
	bool match_on_length)
{

	utility::vector1< bool > cdrs = c;
	if ( cdrs.size() == 6 ) {
		cdrs.push_back(false);
		cdrs.push_back(false);
	}

	CDRDBSequenceSet sequence_set;

	for ( core::Size i = 1; i <= core::Size(ab_info_->get_total_num_CDRs( true /* include CDR4 */)); ++i ) {
		if ( ! cdrs[ i ] ) { continue; }
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );

		if ( cdr == l4 || cdr == h4 ) {
			TR << "Skipping L4/H4 sequence loading.  No profiles exist!" << std::endl;
			continue;
		}
		// Load data from either the stored DataCache or AntibodyInfo.
		CDRClusterEnum current_cluster;

		if ( pose.data().has(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO) ) {
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
			current_cluster = cluster_cache.get_cluster(cdr)->cluster();
		} else {
			current_cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
		}
		if ( current_cluster  == NA ) {
			TR<< ab_info_->get_CDR_name(cdr) << " : Unable to identify cluster for cdr.  Skipping..." << std::endl;
			continue;

		}

		sequence_set[ cdr ] = load_cdr_sequences_for_cdr(
			cdr,
			ab_info_->get_CDR_length(cdr, pose, North),
			current_cluster,
			match_on_length);
	}
	return sequence_set;
}

utility::vector1<CDRDBSequence>
AntibodyDatabaseManager::load_cdr_sequences_for_cdr_length(
	const CDRNameEnum cdr,
	const core::Size length)
{
	//Use a bogus Cluster
	return load_cdr_sequences_for_cdr(cdr, length, static_cast<CDRClusterEnum>(1), true);
}

utility::vector1<CDRDBSequence>
AntibodyDatabaseManager::load_cdr_sequences_for_cdr_cluster(
	const CDRNameEnum cdr,
	const clusters::CDRClusterEnum cluster)
{
	//Use a bogus length
	return load_cdr_sequences_for_cdr(cdr, 1, cluster, false);
}

utility::vector1<CDRDBSequence>
AntibodyDatabaseManager::load_cdr_sequences_for_cdr(
	CDRNameEnum const cdr,
	core::Size const length,
	CDRClusterEnum const cluster,
	bool load_on_length)
{

	utility::vector1<CDRDBSequence> sequence_set;

	std::string base_statement =
		"SELECT\n"
		"\tseq,\n "
		"\tPDB,\n"
		"\tfullcluster,\n"
		"\toriginal_chain,\n"
		"\tdis,\n"
		"\tDistDegree,\n"
		"\tresolution,\n"
		"\trama\n"
		"FROM\n"
		"\tcdr_data \n"
		"WHERE\n"
		"\tdatatag != ? AND"
		"\tCDR = ?";

	//Optional part of Query
	bool use_outliers;
	if ( cdr == h3 ) {
		use_outliers = use_h3_graft_outliers_;
	} else {
		use_outliers = use_outliers_;
	}

	if ( ! use_outliers ) {
		base_statement = base_statement +"\n"+
			"          AND (bb_rmsd_cdr_align < 1.5 OR DistDegree < 40)\n"
			"          AND DistDegree != -1\n"
			"          AND bb_rmsd_cdr_align != -1\n";
	}

	if ( load_on_length ) {
		base_statement = base_statement + "\n"+
			"\tAND length = ?\n";
	} else {
		base_statement = base_statement +"\n" +
			"\tAND fullcluster = ?\n";
	}

	//Lambda vs Kappa
	if ( use_light_chain_type_ && ab_info_->get_light_chain_type_enum() != unknown && (ab_info_->get_light_chain_type_enum() != lambda) ) {
		base_statement += " AND (gene = 'heavy' OR gene = '"+ab_info_->get_light_chain_type()+"')";
	} else if ( use_light_chain_type_ && ab_info_->get_light_chain_type_enum() == lambda ) {
		base_statement += " AND (gene = 'heavy' OR gene = 'lambda6' OR gene = 'lambda')";
	}

	//TR<< "Final Statement:\n"<< base_statement << std::endl;


	cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, db_session_));
	select_statement.bind(1, "loopKeyNotInPaper");
	select_statement.bind(2, ab_info_->get_CDR_name(cdr));

	if ( load_on_length ) {
		select_statement.bind(3, length);
	} else {
		select_statement.bind(3, ab_info_->get_cluster_name(cluster));
	}


	//Get the data.
	cppdb::result pdbid_result(basic::database::safely_read_from_database(select_statement));


	CDRClusterEnumManagerCOP manager = ab_info_->get_cdr_cluster_enum_manager();

	core::Size non_kinked_removed = 0;

	while ( pdbid_result.next() ) {
		std::string seq;
		std::string pdbid;
		std::string cluster;
		std::string original_chain;
		core::Real dis;
		core::Real normalized_dis;
		core::Real resolution;
		std::string rama;

		pdbid_result >> seq >> pdbid >> cluster >> original_chain >> dis >> normalized_dis >> resolution >> rama;
		//TR << seq << " " << pdbid << " "<< cluster << " " << original_chain << " "<<rama << std::endl;

		// Skip if not kinked and we only want kinked here!
		if ( cdr == h3 && use_only_H3_kinked_ && ! is_H3_rama_kinked(rama) ) {
			non_kinked_removed+=1;
			continue;
		}

		//TR << pdbid<<" "<<cluster << original_chain << std::endl;
		//Convert from my database to Rosetta input
		boost::to_lower(pdbid); //Original structures are lower case to distinguish chain easily in filename.

		CDRClusterEnum cluster_enum;
		if ( manager->cdr_cluster_is_present(cluster) ) {

			cluster_enum = ab_info_->get_cluster_enum(cluster);


			//TR.Debug << "PDBId: " << pdbid<<std::endl;
		} else {
			TR<< cluster << " Not Present in enums! Skipping" << std::endl;
			continue;
		}
		std::string pdbid_in = pdbid+original_chain+"_"+ab_info_->get_CDR_name(cdr)+"_0001";

		CDRDBSequence cdr_seq = CDRDBSequence(seq, cluster_enum, pdbid_in, dis, normalized_dis, resolution);
		sequence_set.push_back( cdr_seq );
	}

	if ( non_kinked_removed > 0 ) {
		TR << "Removed "<< non_kinked_removed << " non-kinked H3 structures" << std::endl;
	}

	TR << "Total "<<ab_info_->get_CDR_name(cdr) << " CDRDBPoses: " << sequence_set.size() << std::endl;
	return sequence_set;
}

void
AntibodyDatabaseManager::check_for_graft_instruction_inconsistencies(AntibodyCDRSetOptions const & instructions) {

	//Double check to make sure everything is kosher before grabbing Poses.
	for ( core::Size i=1; i<=CDRNameEnum_total; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRSetOptionsCOP options = instructions[cdr];
		if ( has_vec_inconsistency(options->include_only_clusters(), options->exclude_clusters()) ) {
			utility_exit_with_message("Inconsistency in CLUSTERs option. Cannot leave out and include option.");
		}

		if ( has_vec_inconsistency(options->include_only_pdbs(), options->exclude_pdbs()) ) {
			utility_exit_with_message("Inconsistency in PDBID option. Cannot leave out and include option.");
		}

		if ( has_vec_inconsistency(options->include_only_species(), options->exclude_species()) ) {
			utility_exit_with_message("Inconsistency in SPECIES option. Cannot leave out and include option.");
		}

		if ( has_vec_inconsistency(options->include_only_germlines(), options->exclude_germlines() ) ) {
			utility_exit_with_message("Inconsistency in GERMLINE option. Cannot leave out and include option.");
		}
	}
}


template < typename T >
bool
AntibodyDatabaseManager::has_vec_inconsistency(vector1<T> const &  include, vector1<T> const & exclude) const  {

	if ( (include.size() >= 1) && (exclude.size() >= 1) ) {
		for ( core::Size j=1; j<=include.size(); ++j ) {
			for ( core::Size k=1; k<=exclude.size(); ++k ) {

				if ( include[j]== exclude[k] ) {
					return true;
				}
			}
		}
	}

	return false;
}

void
AntibodyDatabaseManager::bind_vec_constraint(vector1<std::string> const & vec, cppdb::statement& select_statement, core::Size& col) const {

	//Could not get explicit template function specialization to work
	//TR << "start col " << col << std::endl;
	if ( vec.size() >= 1 ) {
		for ( core::Size j = 1; j <= vec.size(); ++j ) {
			col +=1;
			//TR << "col " << col << std::endl;
			//TR << "Binding: " << vec[j] << std::endl;
			select_statement.bind(col, vec[j]);
		}
	}
	//TR << "end col " << col << std::endl;
}

vector1<std::string>
AntibodyDatabaseManager::get_cluster_string_vec(const vector1<CDRClusterEnum>& clusters){
	vector1<std::string> clusters_str;
	for ( core::Size i = 1; i <= clusters.size(); ++i ) {
		clusters_str.push_back(ab_info_->get_cluster_name(clusters[i]));
		//TR << ab_info_->get_cluster_name(clusters[i]) << std::endl;
	}
	return clusters_str;
}

} //antibody
} //protocols
