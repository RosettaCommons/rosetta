// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <utility/sql_database/DatabaseSessionManager.hh>
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
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/rotamer_recovery/RRReporterSQLite.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/Tracer.fwd.hh>


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
			FeaturesReporterOP( new protocols::features::ResidueFeatures() ));
		rr_features->add_features_reporter(
			FeaturesReporterOP( new protocols::features::HBondFeatures(scfxn) ));

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
