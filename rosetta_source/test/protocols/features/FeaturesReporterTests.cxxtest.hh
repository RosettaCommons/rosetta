// -*- mode:c++;tab-width:2;indent-tabs-mode:nil;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RotamerRecovery.cxxtest.hh
/// @brief  Test RotamerRecovery class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/features/AtomAtomPairFeatures.hh>
#include <protocols/features/GeometricSolvationFeatures.hh>
#include <protocols/features/HBondFeatures.hh>
#include <protocols/features/HBondParameterFeatures.hh>
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
#include <protocols/features/ResidueTypesFeatures.hh>
#include <protocols/features/ResidueBurialFeatures.hh>
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>
#include <protocols/features/RotamerBoltzmannWeightFeatures.hh>
#include <protocols/features/RotamerRecoveryFeatures.hh>
// AUTO-REMOVED #include <protocols/features/SaltBridgeFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>


// Project Headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/file_sys_util.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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

		using core::scoring::getScoreFunction;
		using utility::sql_database::DatabaseSessionManager;
		using namespace protocols::features;
		std::string database_filename("features_reporter_tests.db3");

		core_init();
		//Need this to run the features reporter. Adds orbitals to residues
		basic::options::option[ basic::options::OptionKeys::in::add_orbitals](true);
		pose_1ten_ = fullatom_poseop_from_string( pdb_string_1ten() );
		score_function_ = getScoreFunction();

		score_function_->score(*pose_1ten_);

		features_reporters_.push_back(new AtomAtomPairFeatures());
		features_reporters_.push_back(new GeometricSolvationFeatures());
		features_reporters_.push_back(new HBondFeatures(score_function_));
		features_reporters_.push_back(new HBondParameterFeatures(score_function_));
		features_reporters_.push_back(new OrbitalsFeatures());
		features_reporters_.push_back(new PairFeatures());
		features_reporters_.push_back(new PdbDataFeatures());
		features_reporters_.push_back(new PoseCommentsFeatures());
		features_reporters_.push_back(new PoseConformationFeatures());
		features_reporters_.push_back(new ProtocolFeatures());
		features_reporters_.push_back(new ProteinBackboneTorsionAngleFeatures());
		features_reporters_.push_back(new ProteinBackboneAtomAtomPairFeatures());
		features_reporters_.push_back(new ProteinResidueConformationFeatures());
		features_reporters_.push_back(new ProteinRMSDFeatures(pose_1ten_));
		features_reporters_.push_back(new RadiusOfGyrationFeatures());
		features_reporters_.push_back(new ResidueFeatures(score_function_));
		features_reporters_.push_back(new ResidueTypesFeatures());
		features_reporters_.push_back(new ResidueBurialFeatures());
		features_reporters_.push_back(new ResidueSecondaryStructureFeatures());
		features_reporters_.push_back(new RotamerBoltzmannWeightFeatures(score_function_));
		features_reporters_.push_back(new RotamerRecoveryFeatures(score_function_));
		features_reporters_.push_back(new StructureFeatures());
		features_reporters_.push_back(new StructureScoresFeatures(score_function_));

		utility::file::file_delete(database_filename);

		db_session_ = DatabaseSessionManager::get_instance()->get_session(
			"features_reporter_tests.db3");
	}

	void test_RotamerRecovery_main() {
		do_test_schema();

	}

	void do_test_schema() {
		using protocols::features::FeaturesReporterOP;

		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			tr << "Writing schema for '" << reporter->type_name() << "'" << std::endl;
			cppdb::statement schema = (*db_session_) << reporter->schema();
			schema.exec();
		}
	}

	void do_test_report_features() {
		using protocols::features::FeaturesReporterOP;

		core::Size fake_parent_id = 0;
		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			tr << "Reporting features for '" << reporter->type_name() << "'" << std::endl;
			reporter->report_features(*pose_1ten_, fake_parent_id, db_session_);
		}
	}

private:
	core::pose::PoseOP pose_1ten_;
	core::scoring::ScoreFunctionOP score_function_;
	utility::vector1< protocols::features::FeaturesReporterOP > features_reporters_;

	utility::sql_database::sessionOP db_session_;

};
