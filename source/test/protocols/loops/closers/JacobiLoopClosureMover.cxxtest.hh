// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/loops/closers/JacobiLoopClosureMover.cxxtest.hh
/// @brief  Unit tests to verify the JacobiLoopClosureMover and its subfunctions
/// @author teunhoevenaars (teunhoevenaars@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.hh>
#include <protocols/moves/MoverFactory.hh>

// Core Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR("JacobiLoopClosureMoverTests");


class JacobiLoopClosureMoverTests : public CxxTest::TestSuite {
	//Define Variables

private:
	core::pose::PoseOP pose_;

public:

	void setUp(){
		core_init();
		pose_ = core::import_pose::pose_from_file( "protocols/loops/2GB3.pdb" , core::import_pose::PDB_file);
	}

	// Test the default constructor
	void test_default_constructor(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// initialize empty Jacobi Loop Closure Mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover;

		// apply mover without a loop and mm, this should throw an exception
		TS_ASSERT_THROWS(jacobi_mover.apply(pose), const std::exception &);
	}

	// Test the Mover initialisation, and the homogeneous transformations used to define the target and current
	// conformation of the loop
	void test_loop_without_cut(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 7, 14, 14 );

		// initialize Jacobi Loop Closure Mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop);

		// apply mover without any modifications to the pose. The target and current pose should be the same, and therefore
		// the protocol should complete without going into the while-loop of the closure protocol. I.e. finish after 0 cycles
		jacobi_mover.apply(pose);

		// assert that protocol finished successfully after 0 cycles
		TS_ASSERT( jacobi_mover.last_closure_success() == true);
		TS_ASSERT_EQUALS(jacobi_mover.last_closure_cycles(), 0 );
	}

	// Test exception to catch bad loop input
	void test_loop_refinement_inverted(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 14, 7, 10 );

		// initialize Jacobi Loop Closure Mover
		TS_ASSERT_THROWS(protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop), const std::exception & );
	}

	// Test exception to catch bad foldtree in the pose
	void test_loop_with_bad_foldtree(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 7, 14, 10 );

		// initialize Jacobi Loop Closure Mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop);

		// set up FoldTree with a jump point in the loop
		core::kinematics::FoldTree ft = pose.fold_tree();
		ft.new_jump( loop.start()-1, loop.stop()-2, loop.cut());
		ft.reorder( loop.start()-1 );
		pose.fold_tree( ft );

		// apply mover with faulty input pose and catch exception
		TS_ASSERT_THROWS(jacobi_mover.apply(pose), const std::exception & );
	}

	// Test a simple loop closure with few DoFs
	void test_loop_refinement_basic(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 7, 14, 14 );

		// initialize Jacobi Loop Closure Mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop);

		protocols::loops::set_single_loop_fold_tree( pose, loop );
		protocols::loops::add_single_cutpoint_variant( pose, loop );

		pose.set_phi( 8, pose.phi(8) + 5 );
		pose.set_psi( 12, pose.psi(12) + 5 );

		// apply mover
		jacobi_mover.apply(pose);

		// assert that protocol finished successfully after 7 cycles (verified number)
		TS_ASSERT( jacobi_mover.last_closure_success() == true );
		TS_ASSERT_EQUALS (jacobi_mover.last_closure_cycles(), 5);
	}

	// Test a simple loop closure with few DoFs, but cutpoint in the middle of the loop
	void test_loop_refinement_cut_middle(){
		// create variable pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 7, 14, 10 );

		// initialize Jacobi Loop Closure Mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop);

		protocols::loops::set_single_loop_fold_tree( pose, loop );
		protocols::loops::add_single_cutpoint_variant( pose, loop );

		pose.set_phi( 8, pose.phi(8) + 5 );
		pose.set_psi( 12, pose.psi(12) + 5 );

		// apply mover
		jacobi_mover.apply(pose);

		// assert that protocol finished successfully after 10 cycles (verified number)
		TS_ASSERT( jacobi_mover.last_closure_success() == true );
		TS_ASSERT_EQUALS (jacobi_mover.last_closure_cycles(), 7);
	}

	// Test a more advanced loop closure with many DoFs, of which part is fixed via a MoveMap
	void test_loop_reconstruction(){
		// create pose and original pose
		core::pose::Pose pose(*pose_);

		// create sample loop
		protocols::loops::Loop const loop( 7, 37, 14 );

		// create MoveMap
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( true ); // start with all motions allowed
		for ( core::Size i = 1; i < pose.size(); ++i ) { // for all residues
			mm->set( core::id::TorsionID( i, core::id::BB, 3 ), false ); // always fix omega
			if ( (i > 15 && i < 25) || i > 34 ) { // fix phi and psi angles for subset
				mm->set( core::id::TorsionID( i, core::id::BB, 1 ), false );
				mm->set( core::id::TorsionID( i, core::id::BB, 2 ), false );
			}
		}

		// initialize Jacobi Loop Closure Mover for complex move
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover jacobi_mover (loop, mm);
		jacobi_mover.set_dq_max_allowed(numeric::constants::d::pi_over_2);

		protocols::loops::set_single_loop_fold_tree( pose, loop );
		protocols::loops::add_single_cutpoint_variant( pose, loop );

		// disrupt loop. N.B. all torsion angles are reset, including those fixed in the mover using the MoveMap
		for ( core::Size i = loop.start(); i < loop.stop(); ++i ) {
			pose.set_phi( i, 50 );
			pose.set_psi( i, 100 );
		}

		// store pose before mover is applied
		core::pose::Pose original_pose(pose);

		// apply mover
		jacobi_mover.apply(pose);

		// assert that protocol finished successfully after 18 cycles (verified number)
		TS_ASSERT( jacobi_mover.last_closure_success() );
		TS_ASSERT_EQUALS (jacobi_mover.last_closure_cycles(), 27);

		// assert that torsion have changed as expected
		for ( core::Size i = 1; i < pose.size(); ++i ) { // check all residues
			for ( core::Size j = 1; j <= 3; ++j ) { // check torsion angles
				// define current TorsionID
				core::id::TorsionID check_torsionid(i, core::id::BB, j);
				// check that all TorsionIDs that are inside MoveMap and loop are different after move
				if ( mm -> get(check_torsionid) && i >= loop.start() && i <= loop.stop() ) {
					// don't check last two residues in the loop because its torsions are changed when the jump is removed from the Pose
					// absolute values to avoid issues with +180 and -180.
					TS_ASSERT_DIFFERS(std::abs(original_pose.torsion(check_torsionid)),
						std::abs(pose.torsion(check_torsionid)) );
				}
				// check that all TorsionIDs that are outside MoveMap have remained the same
				if ( !mm -> get(check_torsionid) && ! ( i == loop.stop() || i == loop.stop()+1 ) ) {
					// don't check last two residues in the loop because its torsions are changed when the jump is removed from the Pose
					// absolute values to avoid issues with +180 and -180.
					TS_ASSERT_DELTA(std::abs(original_pose.torsion(check_torsionid)),
						std::abs(pose.torsion(check_torsionid)), 1e-3 );
				}
			}
		}
	}

	// test the copy constructor to make sure we are getting deep copies.
	void test_JacobiLoopClosureMover_copy_ctor() {
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover mover1;

		// Set some values in m1 to non-defaults.
		mover1.set_max_cycles( 33 );
		mover1.set_error_norms(3e-3, 3.3e-3);

		// copy mover1
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover mover2( mover1 );

		// Make sure they start off the same
		TS_ASSERT_EQUALS( mover1.get_max_cycles(), mover2.get_max_cycles() );
		TS_ASSERT_EQUALS( mover1.get_allowed_norm_lin(), mover2.get_allowed_norm_lin() );
		TS_ASSERT_EQUALS( mover1.get_allowed_norm_rot(), mover2.get_allowed_norm_rot() );
		TS_ASSERT_EQUALS( mover1.get_dq_max_allowed(), mover2.get_dq_max_allowed() );

		// change the values in m2
		mover2.set_max_cycles( 34 );
		mover2.set_error_norms(4e-3, 3.4e-3);
		mover2.set_dq_max_allowed(1.0);

		// Make sure they are different now
		TS_ASSERT( mover1.get_max_cycles() != mover2.get_max_cycles() );
		TS_ASSERT( mover1.get_allowed_norm_lin() != mover2.get_allowed_norm_lin() );
		TS_ASSERT( mover1.get_allowed_norm_rot() != mover2.get_allowed_norm_rot() );
		TS_ASSERT( mover1.get_dq_max_allowed() != mover2.get_dq_max_allowed() );
	}

	// test the MoverFactory to make sure we can create the mover and down cast appropriately
	void test_JacobiLoopClosureMover_mover_factory() {
		std::string mover_name = "JacobiLoopClosureMover";
		protocols::moves::MoverFactory * mover_factory = protocols::moves::MoverFactory::get_instance();
		protocols::moves::MoverOP my_mover = mover_factory->newMover( mover_name );

		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMoverOP jacobi_mover = utility::pointer::dynamic_pointer_cast< protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover > ( my_mover );
		TS_ASSERT( jacobi_mover ); // make sure we got back the right mover type

		// Create a Tag instance to test for RosettaScripts compatibility
		basic::datacache::DataMap dm;

		using namespace utility::tag;
		using core::Real;

		// Set up tag so all options are not default
		TagOP tag ( utility::pointer::make_shared< Tag >() );
		tag->setName( mover_name );
		tag->setOption< Real >( "max_cycles", 0 );
		tag->setOption< Real >( "err_rot_allowed", 0 );
		tag->setOption< Real >( "err_lin_allowed", 0 );
		tag->setOption< bool >( "verbose", false);
		// set up Loop subtag
		TagOP subtag ( utility::pointer::make_shared< Tag >() );
		subtag->setName( "Loop" );
		subtag->setOption< core::Size >( "start", 10 );
		subtag->setOption< core::Size >( "stop", 13 );
		tag->addTag(subtag);

		// Only the Tag is used by the JacobiLoopClosureMover's parse_my_tag, so the other instances that I'm passing through
		// can be nonsense
		protocols::moves::MoverOP my_configured_mover = mover_factory->newMover( tag, dm );

		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMoverOP configured_jacobi_mover = utility::pointer::dynamic_pointer_cast< protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover > ( my_configured_mover );
		TS_ASSERT( configured_jacobi_mover ); // make sure we got back the right mover type

		// Make sure the mover's configuration reflect what is in the Tag (no function implemented to read verbose)
		TS_ASSERT_EQUALS( configured_jacobi_mover->get_max_cycles(), 0 );
		TS_ASSERT_EQUALS( configured_jacobi_mover->get_allowed_norm_lin(), 0 );
		TS_ASSERT_EQUALS( configured_jacobi_mover->get_allowed_norm_rot(), 0 );
		TS_ASSERT_EQUALS( configured_jacobi_mover->get_loop().start(), 10 );
		TS_ASSERT_EQUALS( configured_jacobi_mover->get_loop().stop(), 13 );
	}

	void tearDown(){

	}
};
