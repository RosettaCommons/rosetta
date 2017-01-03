// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/movers/LoopLengthChange.cxxtest.hh
/// @brief  test for LoopLengthChange mover
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>

#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>

#include <core/pose/annotated_sequence.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("test.protocols.protein_interface_design.movers.LoopLengthChange.cxxtest.hh");

// --------------- Test Class --------------- //

class LoopLengthChangeTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_insert_one_then_three() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "SML", "fa_standard");

		TR << "making mover..." << std::endl;

		protocols::protein_interface_design::movers::LoopLengthChange llcm;
		llcm.loop_start(1);
		llcm.loop_end(2);
		llcm.loop_cut(1);
		llcm.delta(1);

		TR << "applying mover..." << std::endl;
		llcm.apply(pose);

		TR << "testing 1..." << std::endl;
		//did it grow by 1, inserting A?
		TS_ASSERT_EQUALS( pose.size(), 4);
		TS_ASSERT_EQUALS( pose.residue_type(3).name1(), 'A');  //why it inserts after the end of the loop is...a design flaw?

		TR << "altering mover..." << std::endl;
		//change to insert 3 Ws
		llcm.restype_char('W');
		llcm.delta(3);

		TR << "reapplying mover..." << std::endl;
		llcm.apply(pose);

		TR << "testing 2..." << std::endl;
		//did it grow by 3, inserting W?
		TS_ASSERT_EQUALS( pose.size(), 7);
		TS_ASSERT_EQUALS( pose.residue_type(3).name1(), 'W');
		TS_ASSERT_EQUALS( pose.residue_type(4).name1(), 'W');
		TS_ASSERT_EQUALS( pose.residue_type(5).name1(), 'W');
	}

};
