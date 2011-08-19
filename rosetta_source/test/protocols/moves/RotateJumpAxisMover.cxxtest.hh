// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/RotateJumpAxisMover.cxxtest.hh
/// @brief  test for RotateJumpAxisMover
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header
#include <protocols/moves/RotateJumpAxisMover.hh>

// project headers
#include <core/types.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/Stub.hh>
//#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


// --------------- Test Class --------------- //

class RotateJumpAxisMoverTests : public CxxTest::TestSuite {

	core::pose::Pose pose; //RotateJumpAxisMover.pdb
	core::Size rb_jump; //1

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb( pose, "protocols/moves/RotateJumpAxisMover.pdb" );
		rb_jump = 1;

		//we need a pose with a jump in it.  The input pose has 6 residues, we'll tear it into two peices
		//1-3 and 4-6
		using core::kinematics::Edge;
		Edge const jump_edge(2, 5, 1, "CA", "CA", false);
		core::kinematics::FoldTree fold_tree(pose.total_residue());
		fold_tree.add_edge(Edge(1, 2, Edge::PEPTIDE));
		fold_tree.add_edge(Edge(2, 3, Edge::PEPTIDE));
		fold_tree.add_edge(jump_edge);

		fold_tree.add_edge(Edge(4, 5, Edge::PEPTIDE));
		fold_tree.add_edge(Edge(5, pose.total_residue(), Edge::PEPTIDE));

		fold_tree.delete_unordered_edge(1, pose.total_residue(), Edge::PEPTIDE);
		fold_tree.reorder(1);
		pose.fold_tree(fold_tree);
		//std::cout << "set_up() fold tree: " << pose.fold_tree() << std::endl;
	}

	void tearDown() {
		pose.clear();
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	///@brief test a RJAMover running in single angle mode
	void test_RotateJumpAxisMover_single() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::moves::RotateJumpAxisMoverOP;
		using protocols::moves::RotateJumpAxisMover;
		core::Angle degrees(10);
		RotateJumpAxisMoverOP RJA_mover = new RotateJumpAxisMover(rb_jump, degrees);

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/moves/RotateJumpAxisMover_single.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_single

	///@brief test a RJAMover running in range mode
	void test_RotateJumpAxisMover_range() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::moves::RotateJumpAxisMoverOP;
		using protocols::moves::RotateJumpAxisMover;
		core::Angle upper(20), lower(10);
		RotateJumpAxisMoverOP RJA_mover = new RotateJumpAxisMover(rb_jump, lower, upper);

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/moves/RotateJumpAxisMover_range.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_range

	///@brief test a RJAMover running in random mode
	void test_RotateJumpAxisMover_random() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::moves::RotateJumpAxisMoverOP;
		using protocols::moves::RotateJumpAxisMover;
		RotateJumpAxisMoverOP RJA_mover = new RotateJumpAxisMover(rb_jump);

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/moves/RotateJumpAxisMover_random.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_random

};//end class
