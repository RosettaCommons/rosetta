// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/SandwichFeatures.cxxtest.hh
/// @brief  Test SandwichFeatures
/// @author Doonam Kim (doonam.kim@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>
#include <protocols/features/SecondaryStructureSegmentFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/StructureFeatures.hh>

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

#include <utility/tag/Tag.fwd.hh> // for parse_my_tag

static basic::Tracer tr("protocols.features.SandwichFeaturesTests.cxxtest");

class SandwichFeaturesTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
  	tr << "setUp ()" << std::endl;
		tr << "which version I compiled : 	 ss << " << std::endl;
		using utility::sql_database::DatabaseSessionManager;
		using namespace protocols::features;

		core_init();

		database_filename_ = "SandwichFeatures_reporter_tests.db3";
		utility::file::file_delete(database_filename_);
		db_session_OP = basic::database::get_db_session(database_filename_);

		pose_3B83_OP = core::import_pose::pose_from_pdb("protocols/features/3B83_nochain.pdb");
    relevant_residues_3B83_ = utility::vector1< bool >(pose_3B83_OP->total_residue(), true);
		batch_id_ = 0;

		structure_reporter_ = protocols::features::StructureFeaturesOP( new StructureFeatures() );

		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new ResidueSecondaryStructureFeatures() ));
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new SecondaryStructureSegmentFeatures() )); // it is needed for sandwichfeatures
		features_reporters_.push_back(protocols::features::FeaturesReporterOP( new strand_assembly::SandwichFeatures() ));

    // Use arbitrary outputtag and inputtag, instead of being dependent on getting this from JD2
    output_tag_ = "FeaturesReporterTests_outputtag";
    input_tag_ = "FeaturesReporterTests_inputtag";

		write_full_schema(db_session_OP);

	}

	void test_main() {
		report_ResidueSecondaryStructureFeatures();
		report_SecondaryStructureSegmentFeatures();
		report_SandwichFeatures();
	}

	void write_full_schema(utility::sql_database::sessionOP db_session) {
		tr << "write_full_schema" << std::endl;
		using protocols::features::FeaturesReporterOP;

		structure_reporter_->write_schema_to_db(db_session);

		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			tr << "Writing schema for '" << reporter->type_name() << "'" << std::endl;
			reporter->write_schema_to_db(db_session);
		}

	}

  void report_ResidueSecondaryStructureFeatures() {
    		tr << "report_ResidueSecondaryStructureFeatures" << std::endl;

				ResidueSecondaryStructureFeatures_reporters_OP = protocols::features::ResidueSecondaryStructureFeaturesOP( new protocols::features::ResidueSecondaryStructureFeatures() );
				ResidueSecondaryStructureFeatures_reporters_OP->report_features(
																										 *pose_3B83_OP, //	core::pose::Pose const & pose,
																										 relevant_residues_3B83_, //		utility::vector1<bool> const & relevant_residues,
																										 1, //parent_id,// 		StructureID struct_id
																										 db_session_OP //	utility::sql_database::sessionOP db_session
																										 );
	}

  void report_SecondaryStructureSegmentFeatures() {
    		tr << "report_SecondaryStructureSegmentFeatures" << std::endl;
				SecondaryStructureSegmentFeatures_reporters_OP = protocols::features::SecondaryStructureSegmentFeaturesOP( new protocols::features::SecondaryStructureSegmentFeatures() );
				SecondaryStructureSegmentFeatures_reporters_OP->report_features(
																										 *pose_3B83_OP, //	core::pose::Pose const & pose,
																										 relevant_residues_3B83_, //		utility::vector1<bool> const & relevant_residues,
																										 1, //parent_id,// 		StructureID struct_id
																										 db_session_OP //	utility::sql_database::sessionOP db_session
																										 );
	}

  void report_SandwichFeatures() {
    		tr << "report_SandwichFeatures" << std::endl;

				tr << "size of relevant_residues_3B83_ vector: " << relevant_residues_3B83_.size() << std::endl;
				TS_ASSERT( relevant_residues_3B83_.size() == 94 );

				SandwichFeatures_reporters_OP = protocols::features::strand_assembly::SandwichFeaturesOP( new protocols::features::strand_assembly::SandwichFeatures() );

				std::string const string_1 = "protot/ubq_frag.pdb"; //bogus
				std::string const string_2 = "sele";
				std::stringstream ss;
				ss << "<RigidChuM name=cunk mplate=\"" << string_1 << "\" sector=\"" << string_2 << "\" />";
				utility::tag::TagPtr tag( new utility::tag::Tag );
				tag->read( ss );

				// to assign initial values
				SandwichFeatures_reporters_OP->parse_my_tag(
																										tag, data, filters, movers, 
																										*pose_3B83_OP //	core::pose::Pose const & /*pose*/);
																										);

				returned_from_report_features = SandwichFeatures_reporters_OP->report_features(
																										 *pose_3B83_OP, //	core::pose::Pose const & pose,
																										 relevant_residues_3B83_, //		utility::vector1<bool> const & relevant_residues,
																										 1, //parent_id,// 		StructureID struct_id
																										 db_session_OP //	utility::sql_database::sessionOP db_session
																										 );

				tr << "returned_from_report_features: " << returned_from_report_features << std::endl;
				TS_ASSERT( returned_from_report_features == 1 );
  }

private:
  core::pose::PoseOP pose_3B83_OP;

	core::scoring::ScoreFunctionOP score_function_;

	utility::vector1<bool> relevant_residues_3B83_;

	protocols::features::StructureFeaturesOP structure_reporter_;

	protocols::features::ResidueSecondaryStructureFeaturesOP ResidueSecondaryStructureFeatures_reporters_OP;
	protocols::features::SecondaryStructureSegmentFeaturesOP SecondaryStructureSegmentFeatures_reporters_OP;
	protocols::features::strand_assembly::SandwichFeaturesOP SandwichFeatures_reporters_OP;

	utility::vector1< protocols::features::FeaturesReporterOP > features_reporters_;

		// for parse_my_tag
		utility::tag::TagCOP tag;
		basic::datacache::DataMap data;
		protocols::filters::Filters_map filters;
		protocols::moves::Movers_map movers;
		// for parse_my_tag

	std::string database_filename_;
	utility::sql_database::sessionOP db_session_OP;
	core::Size batch_id_;

	core::Size returned_from_report_features;

  std::string output_tag_;
  std::string input_tag_;
};
