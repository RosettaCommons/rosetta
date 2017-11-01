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

#include <protocols/protein_interface_design/movers/AddChainBreak.hh>

#include <core/pose/annotated_sequence.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("test.protocols.protein_interface_design.movers.AddChainBreak.cxxtest.hh");

class AddChainBreakTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_pose_manipulation() {
		core::pose::Pose pose;

		// Sanity check no-op
		core::pose::make_pose_from_sequence(pose, "TESTTESTTEST", "fa_standard");
		TS_ASSERT_EQUALS(pose.num_chains(), 1);
		TS_ASSERT_EQUALS(pose.num_jump(), 0);

		protocols::protein_interface_design::movers::AddChainBreak acb;
		acb.find_automatically(true);
		acb.change_foldtree(true);
		acb.change_conformation(true);
		acb.apply(pose);

		TS_ASSERT_EQUALS(pose.num_chains(), 1);
		TS_ASSERT_EQUALS(pose.num_jump(), 0);

		// Just foldtree
		core::pose::make_pose_from_sequence(pose, "TESTTESTTEST", "fa_standard");
		pose.delete_residue_slow(4);
		TS_ASSERT_EQUALS(pose.num_chains(), 1);
		TS_ASSERT_EQUALS(pose.num_jump(), 0);

		acb.find_automatically(true);
		acb.change_foldtree(true);
		acb.change_conformation(false);
		acb.apply(pose);
		TS_ASSERT_EQUALS(pose.num_chains(), 1);
		TS_ASSERT_EQUALS(pose.num_jump(), 1);

		// Just conformation
		core::pose::make_pose_from_sequence(pose, "TESTTESTTEST", "fa_standard");
		pose.delete_residue_slow(4);

		acb.find_automatically(true);
		acb.change_foldtree(false);
		acb.change_conformation(true);
		acb.apply(pose);
		TS_ASSERT_EQUALS(pose.num_chains(), 2);
		TS_ASSERT_EQUALS(pose.num_jump(), 0);

		// Both foldtree & conformation
		core::pose::make_pose_from_sequence(pose, "TESTTESTTEST", "fa_standard");
		pose.delete_residue_slow(4);

		acb.find_automatically(true);
		acb.change_foldtree(true);
		acb.change_conformation(true);
		acb.apply(pose);
		TS_ASSERT_EQUALS(pose.num_chains(), 2);
		TS_ASSERT_EQUALS(pose.num_jump(), 1);

		// Assert that find_automatically is re-entrant
		acb.apply(pose);
		TS_ASSERT_EQUALS(pose.num_chains(), 2);
		TS_ASSERT_EQUALS(pose.num_jump(), 1);
	}

};
