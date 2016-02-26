// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/SimpleThreadingMoverTests
/// @brief  test for SimpleThreadingMover
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/simple_moves/SimpleThreadingMover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>

// Utility Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.simple_moves.SimpleThreadingMoverTests.cxxtest.hh");

// --------------- Test Class --------------- //

class SimpleThreadingMoverTests : public CxxTest::TestSuite {

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pose::Pose ab_pose_aho; //Full PDB
	protocols::antibody::AntibodyInfoOP ab_info_north_aho;

public:

	void setUp() {
		using namespace protocols::antibody;

		core_init();
		core::import_pose::pose_from_file(ab_pose_aho, "protocols/antibody/1bln_AB_aho.pdb");
		ab_info_north_aho = AntibodyInfoOP(new AntibodyInfo(ab_pose_aho, AHO_Scheme, North));
		
	}

	void tearDown() {
	}

	void test_basic_class() {
		using namespace protocols::antibody;
		using namespace protocols::simple_moves;

		core::Size start = ab_info_north_aho->get_CDR_start(l1, ab_pose_aho);
		core::Size end = ab_info_north_aho->get_CDR_end(l1, ab_pose_aho);

		//L1 length is 16, RSSQSIVHSTGNTYLE.
		std::string seq = "ATATATATATATATAT";
		SimpleThreadingMover threader = SimpleThreadingMover(seq, start);
		TS_ASSERT_THROWS_NOTHING( threader.apply(ab_pose_aho) );
		TS_ASSERT_EQUALS(ab_pose_aho.sequence(start, end), seq)


		//With Options


		std::string seq2 = "ASASASASASASASAS";
		threader.set_pack_rounds(3);
		threader.set_pack_neighbors(true);
		threader.set_neighbor_distance(8);
		threader.set_sequence(seq2, start);
		TS_ASSERT_THROWS_NOTHING( threader.apply( ab_pose_aho ));
		TS_ASSERT_EQUALS(ab_pose_aho.sequence(start, end), seq2);

	}
};

