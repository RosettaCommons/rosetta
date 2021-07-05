// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rotamer_recovery/RRReporterSQLite.cxxtest.hh
/// @brief  Reporter Classes for rotamer recovery
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>

// Project Headers
#include <test/core/init_util.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/file/file_sys_util.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/HBondFeatures.hh>
#include <protocols/features/ReportToDB.hh>
// C++ Headers

//Auto Headers
#include <utility/pointer/owning_ptr.hh>
#include <string>


static basic::Tracer TR("protocols.rotamer_recovery.RRReporterSQLite.cxxtest");

class RRReporterSQLiteTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
		core_init();
	}

	void test_RRReporterSQLite_full(){

		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using core::conformation::Residue;
		using namespace basic::database;
		using namespace protocols::rotamer_recovery;

		std::string database_filename(
			"RRReporterSQLite__test_RRReporterSQLite_full.db3");
		utility::file::file_delete(database_filename);

		utility::sql_database::sessionOP db_session(
			get_db_session(database_filename));

		RRReporterSQLite rs;
		rs.set_struct_id1(1);
		rs.set_output_level( protocols::rotamer_recovery::OL_full );
		rs.write_schema_to_db( db_session );
		rs.db_session( db_session );

		Pose pose ( fullatom_pose_from_string( pdb_string_1ten() ) );
		Residue residue1 ( pose.residue(1) );
		Residue residue2 ( pose.residue(2) );
		Real s(0);
		rs.report_rotamer_recovery( pose, pose, residue1, residue1, s++, true );
		rs.report_rotamer_recovery( pose, pose, residue2, residue2, s, false );

		std::string statement_string = "SELECT * FROM rotamer_recovery;";
		cppdb::statement stmt(safely_prepare_statement(statement_string, db_session));
		cppdb::result res(safely_read_from_database(stmt));
		std::string struct1_name, struct2_name;
		Size chain1, res1, chain2, res2;
		std::string name1, name3, residue_type;
		std::string protocol_name, protocol_params, comparer_name, comparer_params;
		Real score;
		Size recovered;


		res.next();
		res >> struct1_name >> chain1 >> res1;
		res >> struct2_name >> chain2 >> res2;
		res >> name1 >> name3 >> residue_type;
		res >> protocol_name >> protocol_params >> comparer_name >> comparer_params;
		res >> score >> recovered;
		TS_ASSERT_EQUALS(struct1_name, "EMPTY_JOB_use_jd2");
		TS_ASSERT_EQUALS(struct2_name, "EMPTY_JOB_use_jd2");
		TS_ASSERT_EQUALS(name1, "L");
		TS_ASSERT_EQUALS(name3, "LEU");
		TS_ASSERT_EQUALS(residue_type, "LEU:NtermProteinFull");
		TS_ASSERT_EQUALS(chain1, 1);
		TS_ASSERT_EQUALS(res1, 1);
		TS_ASSERT_EQUALS(chain2, 1);
		TS_ASSERT_EQUALS(res2, 1);
		TS_ASSERT_EQUALS(protocol_name, "");
		TS_ASSERT_EQUALS(protocol_params, "");
		TS_ASSERT_EQUALS(comparer_name, "");
		TS_ASSERT_EQUALS(comparer_params, "");
		TS_ASSERT_DELTA(score, 0.0, 0.0001);
		TS_ASSERT_EQUALS(recovered, true);

		res.next();
		res >> struct1_name >> chain1 >> res1;
		res >> struct2_name >> chain2 >> res2;
		res >> name1 >> name3 >> residue_type;
		res >> protocol_name >> protocol_params >> comparer_name >> comparer_params;
		res >> score >> recovered;
		TS_ASSERT_EQUALS(struct1_name, "EMPTY_JOB_use_jd2");
		TS_ASSERT_EQUALS(struct2_name, "EMPTY_JOB_use_jd2");
		TS_ASSERT_EQUALS(name1, "D");
		TS_ASSERT_EQUALS(name3, "ASP");
		TS_ASSERT_EQUALS(residue_type, "ASP");
		TS_ASSERT_EQUALS(chain1, 1);
		TS_ASSERT_EQUALS(res1, 2);
		TS_ASSERT_EQUALS(chain2, 1);
		TS_ASSERT_EQUALS(res2, 2);
		TS_ASSERT_EQUALS(protocol_name, "");
		TS_ASSERT_EQUALS(protocol_params, "");
		TS_ASSERT_EQUALS(comparer_name, "");
		TS_ASSERT_EQUALS(comparer_params, "");
		TS_ASSERT_DELTA(score, 1.0, 0.0001);
		TS_ASSERT_EQUALS(recovered, false);

	}

	void test_RRReporterSQLite_features(){

		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using core::conformation::Residue;
		using namespace basic::database;
		using namespace protocols::rotamer_recovery;
		using protocols::features::FeaturesReporterOP;

		core::scoring::ScoreFunctionOP scfxn(core::scoring::get_score_function());
		Pose pose ( fullatom_pose_from_string( pdb_string_1ten() ) );
		scfxn->score(pose);

		std::string features_db_fname(
			"RRReporterSQLite__test_RRReporterSQLite_features_features.db3");
		utility::file::file_delete(features_db_fname);
		utility::sql_database::sessionOP features_db_session(
			get_db_session(features_db_fname));
		protocols::features::ReportToDBOP rr_features( new protocols::features::ReportToDB(
			features_db_session, "rr_features", "Rotamer Recovery Features",
			true, 2000) );

		rr_features->add_features_reporter(
			utility::pointer::make_shared< protocols::features::ResidueFeatures >());
		rr_features->add_features_reporter(
			utility::pointer::make_shared< protocols::features::HBondFeatures >(scfxn));

		std::string results_db_fname(
			"RRReporterSQLite__test_RRReporterSQLite_features_results.db3");
		utility::file::file_delete(results_db_fname);
		utility::sql_database::sessionOP results_db_session(
			get_db_session(results_db_fname));

		RRReporterSQLite rs;
		rs.set_struct_id1(1);
		rs.set_output_level( protocols::rotamer_recovery::OL_features );
		rs.set_predicted_report_to_db(rr_features);
		rs.write_schema_to_db(results_db_session);
		rs.db_session(results_db_session);


		Residue residue1 ( pose.residue(3) );
		Residue residue2 ( pose.residue(75) );
		Real s(0);
		rs.report_rotamer_recovery( pose, pose, residue1, residue1, s, true );
		s += 1;
		rs.report_rotamer_recovery( pose, pose, residue2, residue2, s, false );

		{
			std::string statement_string = "SELECT * FROM rotamer_recovery;";
			cppdb::statement stmt(safely_prepare_statement(
				statement_string, results_db_session));
			cppdb::result res(safely_read_from_database(stmt));
			protocols::features::StructureID struct_id;
			Size resNum;
			Real divergence;
			Size recovered;

			res.next();
			res >> struct_id >> resNum >> divergence >> recovered;
			TS_ASSERT_EQUALS(resNum, 3);
			TS_ASSERT_DELTA(divergence, 0.0, 0.0001);
			TS_ASSERT_EQUALS(recovered, true);

			res.next();
			res >> struct_id >> resNum >> divergence >> recovered;
			TS_ASSERT_EQUALS(resNum, 75);
			TS_ASSERT_DELTA(divergence, 1.0, 0.0001);
			TS_ASSERT_EQUALS(recovered, false);
		}

		{
			std::string statement_string = "SELECT * FROM rotamer_recovery_predictions;";
			cppdb::statement stmt(safely_prepare_statement(
				statement_string, results_db_session));
			cppdb::result res(safely_read_from_database(stmt));
			protocols::features::StructureID struct_id;
			Size resNum;
			protocols::features::StructureID predicted_struct_id;
			Size predicted_resNum;

			res.next();
			res >> struct_id >> resNum >> predicted_struct_id >> predicted_resNum;
			TS_ASSERT_EQUALS(resNum,3);
			TS_ASSERT_EQUALS(predicted_struct_id, 1);
			TS_ASSERT(predicted_resNum = 1); // THIS IS ALMOST CERTAINLY AN ERROR.

			res.next();
			res >> struct_id >> resNum >> predicted_struct_id >> predicted_resNum;
			TS_ASSERT_EQUALS(resNum, 75);
			TS_ASSERT_EQUALS(predicted_struct_id, 2);
			TS_ASSERT_EQUALS(predicted_resNum, 75);
		}

		{
			std::string statement_string = "SELECT count(*) FROM structures;";
			cppdb::statement stmt(safely_prepare_statement(
				statement_string, features_db_session));
			cppdb::result res(safely_read_from_database(stmt));
			Size n_struct;

			res.next();
			res >> n_struct;
			TS_ASSERT_EQUALS(n_struct, 2);
		}
	}


};
