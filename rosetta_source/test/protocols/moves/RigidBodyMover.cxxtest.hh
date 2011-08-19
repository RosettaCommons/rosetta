// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/RigidBodyMover.cxxtest.hh
/// @brief  test for RigidBodyMover
/// @author Sid Chaudhury

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header
#include <protocols/moves/RigidBodyMover.hh>

// project headers
#include <core/types.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/format.hh>


// --------------- Test Class --------------- //

class RigidBodyMoverTests : public CxxTest::TestSuite {

	core::pose::Pose pose; //RigidBodyMover.pdb
	core::Size rb_jump; //1

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb( pose, "protocols/moves/RigidBodyMover.pdb" );
		rb_jump = 1;

		//setting up the fold tree as is used in docking
		core::kinematics::FoldTree fold_tree;

		core::Size jump_pos1 = 197;
		core::Size jump_pos2 = 282;
		core::Size cutpoint = 245;

    fold_tree.clear();
    fold_tree.add_edge( jump_pos1, jump_pos2, rb_jump );
    fold_tree.add_edge( 1, jump_pos1, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos1, cutpoint, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos2, cutpoint+1, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos2, pose.total_residue(), core::kinematics::Edge::PEPTIDE );
    fold_tree.reorder( 1 );

		pose.fold_tree(fold_tree);
		//std::cout << "set_up() fold tree: " << pose.fold_tree() << std::endl;
	}

	void tearDown() {
		pose.clear();
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	///@brief test a RigidBodyPerturbMover without rot_center calculated at the interface
	void test_RigidBodyPerturbMover_no_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodyPerturbMoverOP;
		using protocols::moves::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbMoverOP RB_mover = new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag, protocols::moves::partner_downstream, false);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodyPerturbMover_no_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_no_int

	///@brief test a RigidBodyPerturbMover with rot_center calculated at the interface
	void test_RigidBodyPerturbMover_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodyPerturbMoverOP;
		using protocols::moves::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		core::scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::DOCK_PATCH ) ;
		(*scorefxn)(pose);  //sets up EnergyGraph for interface calculation

		RigidBodyPerturbMoverOP RB_mover = new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag, protocols::moves::partner_downstream, true);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodyPerturbMover_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_int

	///@brief test a RigidBodyPerturbNoCenterMover
	void test_RigidBodyPerturbNoCenterMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodyPerturbNoCenterMoverOP;
		using protocols::moves::RigidBodyPerturbNoCenterMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbNoCenterMoverOP RB_mover = new RigidBodyPerturbNoCenterMover(rb_jump, rot_mag, trans_mag);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodyPerturbNoCenterMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbNoCenterMover

	///@brief test a RigidBodyRandomizeMover
	void test_RigidBodyRandomizeMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodyRandomizeMoverOP;
		using protocols::moves::RigidBodyRandomizeMover;

		RigidBodyRandomizeMoverOP RB_mover = new RigidBodyRandomizeMover(pose, rb_jump, protocols::moves::partner_downstream);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodyRandomizeMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyRandomizeMover

	///@brief test a RigidBodySpinMover
	void test_RigidBodySpinMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodySpinMoverOP;
		using protocols::moves::RigidBodySpinMover;

		RigidBodySpinMoverOP RB_mover = new RigidBodySpinMover(rb_jump);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodySpinMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodySpinMover

	///@brief test a RigidBodyTransMover
	void test_RigidBodyTransMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::moves::RigidBodyTransMoverOP;
		using protocols::moves::RigidBodyTransMover;

		RigidBodyTransMoverOP RB_mover = new RigidBodyTransMover(pose, rb_jump);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/moves/RigidBodyTransMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyTransMover

	///@brief test a RigidBodyTransMover
	//void test_RigidBodyUniformSphereTransMover() {

	//	////////////////////////RBmover///////////////////////////////////////////////
	//	using protocols::moves::UniformSphereTransMoverOP;
	//	using protocols::moves::UniformSphereTransMover;

	//	core::Real mag(10.0);
	//	UniformSphereTransMoverOP RB_mover = new UniformSphereTransMover(rb_jump, mag);

	//	/////////////////////////run
	//	RB_mover->apply(pose);
	//	test::UTracer UT("protocols/moves/RigidBodyUniformSphereTransMover.pdb");
	//	UT.abs_tolerance(0.003);
	//	UT << std::endl;
	//	pose.dump_pdb(UT);

	//}//end test_RigidBodyUniformSphereTransMover

};//end class
