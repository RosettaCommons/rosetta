// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/clusters/ClusterTests.cxxtest.hh
/// @brief  tests for antibody cluster functions
/// @author Jared Adolf-Bryfogle

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/clusters/CDRClusterFeatures.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/FeaturesReporter.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/sql_utils.hh>
#include <cppdb/frontend.h>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>

using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::constraints;
using namespace protocols::features;

static basic::Tracer TR("protocols.antibody.clusters.ClusterTests");

class ClusterTests : public CxxTest::TestSuite {
	core::pose::Pose ab_pose; //Full PDB
	AntibodyInfoOP ab_info;


public:

	void setUp(){
		using utility::sql_database::DatabaseSessionManager;


		core_init();
		core::import_pose::pose_from_file(ab_pose, "protocols/antibody/2j88.pdb", core::import_pose::PDB_file);

		ab_info = AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North) );
		ab_info->setup_CDR_clusters(ab_pose);



	}

	void tearDown(){
		ab_pose.clear();
	}

	void test_features(){
		using protocols::features::StructureID;
		std::string database_filename = "cdr_cluster_features_test.db3";
		utility::file::file_delete(database_filename);
		utility::sql_database::sessionOP db_session = basic::database::get_db_session(database_filename);

		StructureFeaturesOP structure_reporter( new StructureFeatures() );
		FeaturesReporterOP ab_cluster_reporter( new CDRClusterFeatures() );

		structure_reporter->write_schema_to_db(db_session);
		TS_ASSERT_THROWS_NOTHING(ab_cluster_reporter->write_schema_to_db(db_session));

		StructureID parent_id = structure_reporter->report_features(0, db_session, "output_tag", "input_tag");
		TS_ASSERT_THROWS_NOTHING(
			ab_cluster_reporter->report_features(ab_pose, parent_id, db_session)
		);

	}
	void test_cluster_identification(){

		using namespace core::pose::datacache;
		TS_ASSERT(check_if_pose_renumbered_for_clusters(ab_pose));

		TS_ASSERT_EQUALS("L1-11-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l1)->cluster()));
		TS_ASSERT_EQUALS("L2-8-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l2)->cluster()));
		TS_ASSERT_EQUALS("L3-8-1", ab_info->get_cluster_name(ab_info->get_CDR_cluster(l3)->cluster()));

		CDRClusterCOP cluster_L1 = ab_info->get_CDR_cluster(l1);
		TS_ASSERT_EQUALS(l1, cluster_L1->cdr());
		TS_ASSERT_EQUALS('L', cluster_L1->chain());
		TS_ASSERT_EQUALS(24, cluster_L1->pdb_start());
		TS_ASSERT_EQUALS(42, cluster_L1->pdb_end());

		//Need to add numerical values - not sure how to match with possible rounding errors here...
		TS_ASSERT_THROWS_NOTHING(cluster_L1->distance());
		TS_ASSERT_THROWS_NOTHING(cluster_L1->length_normalized_distance());
		TS_ASSERT_THROWS_NOTHING(cluster_L1->normalized_distance_in_degrees());


		TS_ASSERT_EQUALS(l1, ab_info->get_cdr_enum_for_cluster(L1_11_1));
		TS_ASSERT_EQUALS(11, ab_info->get_cluster_length(L1_11_1));
		TS_ASSERT_EQUALS(L1_11_1, ab_info->get_cluster_enum("L1-11-1"));

		ab_info->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(ab_pose);
		TS_ASSERT(ab_pose.data().has(CacheableDataType::CDR_CLUSTER_INFO));
		BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(ab_pose.data().get(CacheableDataType::CDR_CLUSTER_INFO));
		CDRClusterCOP l1_result = cluster_cache.get_cluster(l1);
		CDRClusterCOP l2_result = cluster_cache.get_cluster(l2);
		CDRClusterCOP l3_result = cluster_cache.get_cluster(l3);

		TS_ASSERT_EQUALS(L1_11_1, l1_result->cluster());
		TS_ASSERT_EQUALS(L2_8_1, l2_result->cluster());
		TS_ASSERT_EQUALS(L3_8_1, l3_result->cluster());

		//Test that we can clone the ClusterOP
		CDRClusterCOP l1_cop = ab_info->get_CDR_cluster(l1);
		TS_ASSERT_THROWS_NOTHING(l1_cop->clone());

		//Test manually updating the cluster
		// make sure this is reflected in the datacache after update
		CDRClusterOP new_l1_cluster( new CDRCluster(
			ab_pose,
			l1,
			11,
			L1_11_2,
			ab_info->get_CDR_start(l1, ab_pose),
			.2) );

		ab_info->set_CDR_cluster(l1, new_l1_cluster);
		ab_info->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(ab_pose);
		TS_ASSERT(ab_pose.data().has(CacheableDataType::CDR_CLUSTER_INFO));

		BasicCDRClusterSet const & cluster_cache2 = static_cast< BasicCDRClusterSet const & >(ab_pose.data().get(CacheableDataType::CDR_CLUSTER_INFO));
		CDRClusterCOP new_l1 = cluster_cache2.get_cluster(l1);

		TS_ASSERT_EQUALS(L1_11_2, new_l1->cluster());

	}

};
