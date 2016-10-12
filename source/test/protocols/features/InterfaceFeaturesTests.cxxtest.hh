// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file test/protocols/features/InterfaceFeaturesTests.cxxtest.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/InterfaceFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>

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

static basic::Tracer TR("protocols.features.InterfaceFeaturesTests.cxxtest");

using namespace protocols::features;
using utility::vector1;

class InterfaceFeaturesTests : public CxxTest::TestSuite {

	core::pose::Pose multimer_;
	InterfaceFeaturesOP reporter_;
	std::string db_name_;
	utility::sql_database::sessionOP db_session_;

public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(multimer_, "protocols/features/2J88.pdb", core::import_pose::PDB_file);
		db_name_ =  "InterfaceFeaturesTest.db3";
		reporter_ = InterfaceFeaturesOP( new protocols::features::InterfaceFeatures() );
		utility::file::file_delete(db_name_);
		db_session_ = basic::database::get_db_session(db_name_);
		TR <<"Setup"<<std::endl;
	}

	void tearDown(){
		multimer_.clear();
	}

	void test_reporter(){
		interface_test();
		TR << "Testing features reporter" << std::endl;
		utility::vector1<bool> relavant_residues(multimer_.size(), true);

		reporter_->set_dSASA_cutoff(150); //Not real value
		reporter_->set_pack_separated(false); //speed
		reporter_->set_pack_together(false);

		vector1<std::string> interfaces;
		interfaces.push_back("L_H");
		//interfaces.push_back("L_A");
		//interfaces.push_back("H_A");
		interfaces.push_back("LH_A");

		reporter_->set_interface_chains(interfaces);
		reporter_->set_dSASA_cutoff(150);
		reporter_->set_pack_separated(false);
		reporter_->set_pack_together(false); //speed
		reporter_->set_compute_packstat(false); //speed

		StructureFeaturesOP structure_reporter( new StructureFeatures() );
		structure_reporter->write_schema_to_db(db_session_);
		StructureID parent_id = structure_reporter->report_features(0, db_session_, "output_tag", "input_tag");

		ResidueFeaturesOP residue_reporter( new ResidueFeatures() );
		residue_reporter->write_schema_to_db(db_session_);
		residue_reporter->report_features(multimer_, relavant_residues, parent_id, db_session_);

		TS_ASSERT_THROWS_NOTHING(reporter_->write_schema_to_db(db_session_));
		TS_ASSERT_THROWS_NOTHING(reporter_->report_features(multimer_, relavant_residues, parent_id, db_session_));

	}

	void interface_test(){
		TR << "Testing Interface combos" << std::endl;
		utility::vector1<std::string> interfaces;
		TS_ASSERT_THROWS_NOTHING(reporter_->make_interface_combos(multimer_, interfaces));
		TR << "Interfaces: "<< interfaces.size() << std::endl;
		for ( core::Size i = 1; i<=interfaces.size(); ++i ) {
			TR<<"Interface: " << interfaces[i]<<std::endl;
		}
		TS_ASSERT_EQUALS(interfaces.size(), 6);
	}
};
