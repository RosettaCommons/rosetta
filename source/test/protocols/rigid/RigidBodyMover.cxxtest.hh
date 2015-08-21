// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rigid/RigidBodyMover.cxxtest.hh
/// @brief  test for RigidBodyMover
/// @author Sid Chaudhury

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/rosettascripts.hh>

// Unit header
#include <protocols/rigid/RigidBodyMover.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>


// --------------- Test Class --------------- //

class RigidBodyMoverTests : public CxxTest::TestSuite {

	core::pose::Pose pose; //RigidBodyMover.pdb
	core::Size rb_jump; //1

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb( pose, "protocols/rigid/RigidBodyMover.pdb" );
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
	}

	void tearDown() {
		pose.clear();
	}

	// --------------- Test Cases --------------- //

	/// @brief test a RigidBodyPerturbMover without rot_center calculated at the interface
	void test_RigidBodyPerturbMover_no_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbMoverOP;
		using protocols::rigid::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbMoverOP RB_mover( new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag,
			protocols::rigid::partner_downstream, false) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbMover_no_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_no_int

	/// @brief test a RigidBodyPerturbMover with rot_center calculated at the interface
	void test_RigidBodyPerturbMover_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbMoverOP;
		using protocols::rigid::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		core::scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(
			core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::DOCK_PATCH ) ;
		(*scorefxn)(pose);  //sets up EnergyGraph for interface calculation

		RigidBodyPerturbMoverOP RB_mover( new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag,
			protocols::rigid::partner_downstream, true) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbMover_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_int

	/// @brief test a RigidBodyPerturbNoCenterMover
	void test_RigidBodyPerturbNoCenterMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbNoCenterMoverOP;
		using protocols::rigid::RigidBodyPerturbNoCenterMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbNoCenterMoverOP RB_mover( new RigidBodyPerturbNoCenterMover(rb_jump, rot_mag, trans_mag) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbNoCenterMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbNoCenterMover

	/// @brief test a RigidBodyRandomizeMover
	void test_RigidBodyRandomizeMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyRandomizeMoverOP;
		using protocols::rigid::RigidBodyRandomizeMover;

		RigidBodyRandomizeMoverOP RB_mover( new RigidBodyRandomizeMover(pose, rb_jump,
			protocols::rigid::partner_downstream) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyRandomizeMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyRandomizeMover

	/// @brief test a RigidBodySpinMover
	void test_RigidBodySpinMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodySpinMoverOP;
		using protocols::rigid::RigidBodySpinMover;

		RigidBodySpinMoverOP RB_mover( new RigidBodySpinMover(rb_jump) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodySpinMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodySpinMover

	/// @brief test a RigidBodyTransMover
	void test_RigidBodyTransMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyTransMoverOP;
		using protocols::rigid::RigidBodyTransMover;

		RigidBodyTransMoverOP RB_mover( new RigidBodyTransMover(pose, rb_jump, false) );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyTransMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyTransMover

	/// @brief test parsing RigidBodyTransMover
	void test_RigidBodyTransMover_parse() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyTransMoverOP;
		using protocols::rigid::RigidBodyTransMover;

		RigidBodyTransMoverOP RB_mover( new RigidBodyTransMover() );

		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;
		TagCOP tag = tagptr_from_string("<RigidBodyTransMover name=trans distance=1 jump=1/>\n");
		RB_mover->parse_my_tag( tag, data, filters, movers, pose );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyTransMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyTransMover_parse()
};//end class
