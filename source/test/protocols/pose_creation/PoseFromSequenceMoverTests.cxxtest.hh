// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_creation/PoseFromSequenceMoverTests.cxxtest.hh
/// @brief  Create a pose from a sequence or fasta
/// @author Dan Farrell (danpf@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <basic/datacache/DataMap.hh>
#include <core/types.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <protocols/pose_creation/PoseFromSequenceMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("PoseFromSequenceMoverTests");


class PoseFromSequenceMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:
	PoseFromSequenceMoverTests() {}

	void setUp(){
		core_init();
	}

	void tearDown() {
	}

	void test_ctor() {
		protocols::pose_creation::PoseFromSequenceMover pfsm;
		TS_ASSERT_EQUALS(pfsm.get_residue_type_set(), "fa_standard");
		TS_ASSERT_EQUALS(pfsm.get_extended(), false);
		TS_ASSERT_EQUALS(pfsm.get_sequence(), "");
	}

	void test_parse_my_tag() {
		basic::datacache::DataMap data;
		core::pose::Pose pose(create_test_in_pdb_pose());
		{ // test fasta
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" fasta=\"protocols/pose_creation/PFSM_01.fasta\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_sequence(), "EQLLKALEFLLKELLEKL");
			TS_ASSERT_EQUALS(pfsm.get_residue_type_set(), "fa_standard");
		}
		{ // test fasta with multiple
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" fasta=\"protocols/pose_creation/PFSM_02.fasta\" use_all_in_fasta=\"true\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_sequence(), "EQLLKALEFLLKELLEKL/EQLLKALEFLYYYYYYY/EQLPLERLPRL");
		}
		{ // test fasta with multiple but use only single
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" fasta=\"protocols/pose_creation/PFSM_02.fasta\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_sequence(), "EQLLKALEFLLKELLEKL");
		}
		{ // test sequence  + extended
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" sequence=\"EQLLKALEFLLKELLEKL\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_sequence(), "EQLLKALEFLLKELLEKL");
			TS_ASSERT_EQUALS(pfsm.get_extended(), false);
		}
		{ // test sequence  + extended
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" sequence=\"EQLLKALEFLLKELLEKL\" extended=\"true\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_sequence(), "EQLLKALEFLLKELLEKL");
			TS_ASSERT_EQUALS(pfsm.get_extended(), true);
		}
		{ // test throwing if no seq or fasta
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" extended=\"true\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			TS_ASSERT_THROWS(pfsm.parse_my_tag(tag, data), utility::excn::RosettaScriptsOptionError &);
		}
		{ // test throwing if both sequence and fasta
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" sequence=\"SEVER\" fasta=\"PLYWALL\" extended=\"true\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			TS_ASSERT_THROWS(pfsm.parse_my_tag(tag, data), utility::excn::RosettaScriptsOptionError &);
		}
		{ // test residue type set
			utility::tag::TagCOP tag = tagptr_from_string("<PoseFromSequenceMover name=\"rhino\" sequence=\"SEVER\" residue_type_set=\"fa_standard\" />\n");
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.parse_my_tag(tag, data);
			TS_ASSERT_EQUALS(pfsm.get_residue_type_set(), "fa_standard");
		}
	}

	void test_apply() {
		core::pose::Pose pose;
		{ // don't fail if no seq
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.set_sequence("");
			pfsm.apply(pose);
			TS_ASSERT_EQUALS(pose.size(), 0);
		}
		{ //  basic
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.set_sequence("SEVERPLYWALL");
			pfsm.apply(pose);
			TS_ASSERT_EQUALS(pose.sequence(), "SEVERPLYWALL");
			TS_ASSERT_EQUALS(pose.is_centroid(), false);
		}
		{ //  basic fa
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.set_sequence("SEVERPLYWALL");
			pfsm.set_residue_type_set("centroid");
			pfsm.apply(pose);
			TS_ASSERT_EQUALS(pose.sequence(), "SEVERPLYWALL");
			TS_ASSERT_EQUALS(pose.is_centroid(), true);
		}
		{ //  extended check
			core::pose::Pose extended_pose;
			protocols::pose_creation::PoseFromSequenceMover pfsm;
			pfsm.set_sequence("SEVERPLYWALL");
			pfsm.apply(pose);
			pfsm.set_extended(true);
			pfsm.apply(extended_pose);
			numeric::xyzVector< core::Real > extended_center(core::pose::get_center_of_mass( extended_pose ));
			numeric::xyzVector< core::Real > unextended_center(core::pose::get_center_of_mass( pose ));
			// this is a bit silly but basically just check that the center of mass changes
			TS_ASSERT(extended_center[0] > unextended_center[0] + 5);
			TS_ASSERT(extended_center[1] > unextended_center[1] + 5);
		}
	}
};
