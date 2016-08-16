// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rigid/RotateJumpAxisMover.cxxtest.hh
/// @brief  test for RotateJumpAxisMover
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/pose_funcs.hh>

// Unit header
#include <protocols/rigid/RotateJumpAxisMover.hh>

// Project headers
#include <core/types.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// --------------- Test Class --------------- //

class RotateJumpAxisMoverTests : public CxxTest::TestSuite {

	core::pose::Pose pose;
	core::pose::Pose master_pose; //create_pdb_string_2res_1ten_2res_trp_cage_pose();
	core::Size rb_jump; //1

public:

	// --------------- Fixtures --------------- //
	RotateJumpAxisMoverTests(){
		core_init_with_additional_options("-out:file:no_chainend_ter");
		pose = create_pdb_string_2res_1ten_2res_trp_cage_pose();
		rb_jump = 1;

		pose.dump_pdb("unmodified.pdb");

		//we need a pose with a jump in it.  The input pose has 4 residues, we'll tear it into two peices
		//1-3 and 4-6
		using core::kinematics::Edge;
		Edge const jump_edge(
			1, //residue 1
			3, //to residue 3
			1, //jump labeled 1
			"CA", //atom CA on residue 1
			"CA", //atom CA on residue 3
			false); //IDK
		core::kinematics::FoldTree fold_tree(pose.total_residue());
		fold_tree.clear();
		fold_tree.add_edge(Edge(1, 2, Edge::PEPTIDE));
		fold_tree.add_edge(jump_edge);
		fold_tree.add_edge(Edge(3, pose.total_residue(), Edge::PEPTIDE));
		fold_tree.reorder(1);
		pose.fold_tree(fold_tree);
		//std::cout << "set_up() fold tree: " << pose.fold_tree() << std::endl;

		master_pose = pose;
	}


	virtual ~RotateJumpAxisMoverTests() {}

	static RotateJumpAxisMoverTests *createSuite() {
		return new RotateJumpAxisMoverTests();
	}

	static void destroySuite( RotateJumpAxisMoverTests *suite ) {
		delete suite;
	}

	void setUp() { pose = master_pose; }

	void tearDown() { }


	// --------------- Test Cases --------------- //

	/// @brief test a RJAMover running in single angle mode
	void test_RotateJumpAxisMover_single() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::rigid::RotateJumpAxisMoverOP;
		using protocols::rigid::RotateJumpAxisMover;
		core::Angle degrees(10);
		RotateJumpAxisMoverOP RJA_mover( new RotateJumpAxisMover(rb_jump, degrees) );

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RotateJumpAxisMover_single.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_single

	/// @brief test a RJAMover running in range mode
	void test_RotateJumpAxisMover_range() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::rigid::RotateJumpAxisMoverOP;
		using protocols::rigid::RotateJumpAxisMover;
		core::Angle upper(30), lower(20);
		RotateJumpAxisMoverOP RJA_mover( new RotateJumpAxisMover(rb_jump, lower, upper) );

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RotateJumpAxisMover_range.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_range

	/// @brief test a RJAMover running in random mode
	void test_RotateJumpAxisMover_random() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::rigid::RotateJumpAxisMoverOP;
		using protocols::rigid::RotateJumpAxisMover;
		RotateJumpAxisMoverOP RJA_mover( new RotateJumpAxisMover(rb_jump) );

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RotateJumpAxisMover_random.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_random

	/// @brief test a RJAMover with an underspecifed (not atom-to-atom) Jump
	void test_RotateJumpAxisMover_underspecified_jump() {

		////////////////////////RJAmover///////////////////////////////////////////////
		using protocols::rigid::RotateJumpAxisMoverOP;
		using protocols::rigid::RotateJumpAxisMover;
		core::Angle degrees(10);
		RotateJumpAxisMoverOP RJA_mover( new RotateJumpAxisMover(rb_jump, degrees) );

		//notice this should still produce the same structure as the 10 degree test above!

		//tweak the FoldTree
		core::kinematics::FoldTree tree(pose.total_residue());
		tree.new_jump(1,3,2); //Jump from 1 to 3, cutpoint between 2 and 3.  From a residue level this is the same as the fold tree we already had.
		pose.fold_tree(tree);

		/////////////////////////run
		RJA_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RotateJumpAxisMover_underspecified_jump.pdb");
		pose.dump_pdb(UT);

	}//end test_RotateJumpAxisMover_single

};//end class
