// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RotamerRecovery.cxxtest.hh
/// @brief  Test FeaturesReporter class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/features/AtomAtomPairFeatures.hh>
#include <protocols/features/AtomInResidueAtomInResiduePairFeatures.hh>
#include <protocols/features/AtomTypesFeatures.hh>
#include <protocols/features/BetaTurnDetectionFeatures.hh>
#include <protocols/features/ChargeChargeFeatures.hh>
#include <protocols/features/GeometricSolvationFeatures.hh>
#include <protocols/features/HBondFeatures.hh>
#include <protocols/features/HBondParameterFeatures.hh>
#include <protocols/features/JobDataFeatures.hh>
#include <protocols/features/LoopAnchorFeatures.hh>
#include <protocols/features/OrbitalsFeatures.hh>
#include <protocols/features/PairFeatures.hh>
#include <protocols/features/PdbDataFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinBackboneTorsionAngleFeatures.hh>
#include <protocols/features/ProteinBackboneAtomAtomPairFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ProteinRMSDFeatures.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/RadiusOfGyrationFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/ResidueScoresFeatures.hh>
#include <protocols/features/ResidueTypesFeatures.hh>
#include <protocols/features/ResidueBurialFeatures.hh>
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>
#include <protocols/features/RotamerBoltzmannWeightFeatures.hh>
#include <protocols/features/RotamerRecoveryFeatures.hh>
#include <protocols/features/SaltBridgeFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh> // alphabetically reordered- doonam
#include <protocols/features/strand_assembly/StrandBundleFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>
#include <protocols/features/TrajectoryReportToDB.hh> // alphabetically reordered- doonam
#include <protocols/features/UnrecognizedAtomFeatures.hh>



// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>
#include <utility/file/file_sys_util.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <limits>

// External Headers
#include <cppdb/frontend.h>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>


static basic::Tracer tr("protocols.features.FeaturesReporterTests.cxxtest");

class FeaturesReporterTests : public CxxTest::TestSuite {

public:

	void
	setUp() {

		using core::scoring::get_score_function;
		using utility::sql_database::DatabaseSessionManager;
		using namespace protocols::features;

		core_init();

		database_filename_ = "features_reporter_tests.db3";
		utility::file::file_delete(database_filename_);
		db_session_ = basic::database::get_db_session(database_filename_);

		//Need this to run the features reporter. Adds orbitals to residues
		//basic::options::option[ basic::options::OptionKeys::in::add_orbitals](true); // apl disabling as this screws up the singleton FA_STANDARD residue type set
		//pose_1ten_ = fullatom_poseop_from_string( pdb_string_1ten() );
		pose_1ten_ = core::import_pose::pose_from_pdb("protocols/features/2J88.pdb");
		relevant_residues_ = utility::vector1< bool >(pose_1ten_->total_residue(), true);
		batch_id_ = 0;

		score_function_ = get_score_function();

		score_function_->score(*pose_1ten_);

		structure_reporter_ = protocols::features::StructureFeaturesOP( new StructureFeatures() );

		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new AtomAtomPairFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new AtomInResidueAtomInResiduePairFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new AtomTypesFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new BetaTurnDetectionFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ChargeChargeFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new GeometricSolvationFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new HBondFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new HBondParameterFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new JobDataFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new LoopAnchorFeatures() ));
		//features_reporters_.push_back(new OrbitalsFeatures());
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new PairFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new PdbDataFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new PoseCommentsFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new PoseConformationFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ProteinBackboneTorsionAngleFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ProteinBackboneAtomAtomPairFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ProteinResidueConformationFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ProteinRMSDFeatures(pose_1ten_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new RadiusOfGyrationFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueScoresFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueTypesFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueBurialFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueSecondaryStructureFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new RotamerBoltzmannWeightFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new RotamerRecoveryFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new SaltBridgeFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new strand_assembly::SandwichFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new strand_assembly::StrandBundleFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new StructureScoresFeatures(score_function_) ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new UnrecognizedAtomFeatures() ));

    // Use arbitrary outputtag and inputtag, instead of being dependent on getting this from JD2
    output_tag_ = "FeaturesReporterTests_outputtag";
    input_tag_ = "FeaturesReporterTests_inputtag";
	}

	void test_main() {
		write_full_schema(db_session_);
		do_test_type_name();
		//do_test_report_features(); JAB - tried to enable this so each features reporter is actually tested, but I get an exceeded the timeout error from run.py
	}

	void write_full_schema(utility::sql_database::sessionOP db_session) {
		using protocols::features::FeaturesReporterOP;

		structure_reporter_->write_schema_to_db(db_session);

		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			tr << "Writing schema for '" << reporter->type_name() << "'" << std::endl;
			reporter->write_schema_to_db(db_session);
		}
	}

	void do_test_type_name() {
		using protocols::features::FeaturesReporterOP;
		foreach(FeaturesReporterOP const & reporter, features_reporters_){
			TS_ASSERT_DIFFERS(reporter->type_name(), "Unknown_FeaturesReporter");
		}
	}

	void do_test_report_features() {
		using protocols::features::FeaturesReporterOP;
		using protocols::features::StructureID;

		tr << "Creating structure entry." << std::endl;
		StructureID parent_id = structure_reporter_->report_features(batch_id_, db_session_, output_tag_, input_tag_);
		tr << "Created structure id:" << parent_id << std::endl;

		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			tr << "Reporting features for '" << reporter->type_name() << "'" << std::endl;
			reporter->report_features(*pose_1ten_, batch_id_, db_session_);
		}
	}

	void test_identifier_generation() {
		structure_reporter_->write_schema_to_db(db_session_);
		using protocols::features::FeaturesReporterOP;
		using protocols::features::StructureID;

		tr << "Creating structure entry." << std::endl;

		StructureID struct_id_one = structure_reporter_->report_features(batch_id_, db_session_, output_tag_, input_tag_);
		tr << "Created structure id:" << struct_id_one << std::endl;

		StructureID struct_id_two = structure_reporter_->report_features(batch_id_, db_session_, output_tag_, input_tag_);
		tr << "Created structure id:" << struct_id_two << std::endl;

		TS_ASSERT(struct_id_one == 1);
		TS_ASSERT(struct_id_two == 2);

		// Create a second database session in a separate partition
		// Delete the partitioned sqlite database
		utility::file::file_delete(database_filename_ + "_1");
		utility::sql_database::sessionOP db_session_p1 = utility::sql_database::DatabaseSessionManager::get_instance()->get_session_sqlite3(
			database_filename_,
			utility::sql_database::TransactionMode::standard,
			0,
			false,
			1);

		structure_reporter_->write_schema_to_db(db_session_p1);

		StructureID structure_prefix = 1;
		structure_prefix = structure_prefix << 32;

		tr << "Creating partitioned structure entry in database partition 1, structure prefix: " << structure_prefix << std::endl;

		StructureID struct_id_one_partition1 = structure_reporter_->report_features(batch_id_, db_session_p1, output_tag_, input_tag_);
		tr << "Created structure id:" << struct_id_one_partition1 << std::endl;

		StructureID struct_id_two_partition1 = structure_reporter_->report_features(batch_id_, db_session_p1, output_tag_, input_tag_);
		tr << "Created structure id:" << struct_id_two_partition1 << std::endl;


		TS_ASSERT(struct_id_one_partition1 == 1 + structure_prefix);
		TS_ASSERT(struct_id_two_partition1 == 2 + structure_prefix);
	}

  void test_trajectory_report_to_db() {
    using protocols::features::TrajectoryReportToDBOP;
    TrajectoryReportToDBOP traj_reporter( new protocols::features::TrajectoryReportToDB(
      db_session_, "fake_batch_name", "fake_batch_description",
      false, // false = don't use transactions
      1 // cache_size
    ) );

    // Verify that default stride is 1 (to help make sure changes to this default are noticed)
    TS_ASSERT(traj_reporter->get_stride() == 1);

    TS_ASSERT_THROWS_NOTHING( traj_reporter->set_stride(2) );
    TS_ASSERT( traj_reporter->get_stride() == 2 );

    TS_ASSERT_THROWS_NOTHING( traj_reporter->apply(*pose_1ten_) );

    // Verify that cycle_counts is initialized correctly
    TS_ASSERT( traj_reporter->get_cycle_counts().size() == 1 );
    TS_ASSERT( traj_reporter->get_cycle_counts().begin()->second == 1 );

    // Verify that cycle counting is working correctly for same pose apply
    TS_ASSERT_THROWS_NOTHING( traj_reporter->apply(*pose_1ten_) );
    TS_ASSERT( traj_reporter->get_cycle_counts().begin()->second == 2 );
  }

private:
	core::pose::PoseOP pose_1ten_;
	core::scoring::ScoreFunctionOP score_function_;

	utility::vector1<bool> relevant_residues_;
	core::Size batch_id_;

	protocols::features::StructureFeaturesOP structure_reporter_;
	utility::vector1< protocols::features::FeaturesReporterOP > features_reporters_;

	std::string database_filename_;
	utility::sql_database::sessionOP db_session_;

  std::string output_tag_;
  std::string input_tag_;
};
