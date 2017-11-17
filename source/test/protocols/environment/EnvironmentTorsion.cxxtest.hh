// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/environment/DofPassport.hh>

#include <test/protocols/environment/TorsionMover.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/Environment.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

#include <basic/Tracer.hh>

//C++ headers
#include <iostream>

static basic::Tracer TR("protocols.environment.EnvironmentTorsion.cxxtest");

// --------------- Test Class --------------- //

class EnvironmentTorsion : public CxxTest::TestSuite {
public:

	// Shared data elements go here.
	core::pose::Pose pose;
	utility::vector1< core::Real > init_phis;

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();

		using namespace protocols::environment;
		using namespace core::environment;

		core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

		//store phi values for comparison in the test
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			init_phis.push_back( pose.phi( i ) );
		}

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_empty_environment(){
		TR <<  "Beginning: test_empty_environment"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		TorsionMoverOP mover( new TorsionMover( false, false ) );
		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;

		core::pose::Pose brokered;
		TS_ASSERT_THROWS_NOTHING( brokered = env.start( pose ) );
		core::pose::Pose done;
		TS_ASSERT_THROWS_NOTHING( done = env.end( brokered ) );

		TR <<  "End: test_empty_environment"  << std::endl;
	}

	void test_dual_environment(){
		TR <<  "Beginning: test_dual_environment"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		TorsionMoverOP tier1_mover( new TorsionMover( true, true ) );
		TorsionMoverOP tier2_mover( new TorsionMover( true, true ) );

		EnvironmentOP env1_op( new Environment( "tier1" ) );
		Environment & env1 = *env1_op;
		EnvironmentOP env2_op( new Environment( "tier2" ) );
		Environment & env2 = *env2_op;

		env1.register_mover( tier1_mover );
		env2.register_mover( tier2_mover );

		core::pose::Pose final_pose;

		{
			core::pose::Pose protected1;
			TS_ASSERT_THROWS_NOTHING( protected1 = env1.start( pose ) );
			TS_ASSERT_THROWS( tier2_mover->apply( protected1 ), utility::excn::NullPointerError );
			TS_ASSERT_EQUALS( protected1.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ] );
			{
				core::pose::Pose protected2;
				TS_ASSERT_THROWS_NOTHING( protected2 = env2.start( protected1 ) );
				TS_ASSERT_THROWS( tier1_mover->apply( protected2 ), EXCN_Env_Security_Exception );
				TS_ASSERT_DELTA( protected2.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 0.001 );
				TS_ASSERT_THROWS_NOTHING( tier2_mover->apply( protected2 ) );
				TS_ASSERT_DELTA( protected2.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );
				TS_ASSERT_THROWS_NOTHING( protected1 = env2.end( protected2 ) );
			}
			TS_ASSERT_DELTA( protected1.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );
			TS_ASSERT_THROWS_NOTHING( tier1_mover->apply( protected1) );
			TS_ASSERT_THROWS( tier2_mover->apply( protected1 ), utility::excn::NullPointerError );
			TS_ASSERT_THROWS_NOTHING( final_pose = env1.end( protected1 ) );
		}

		for ( core::Size seqpos = 1; seqpos <= pose.size(); ++seqpos ) {
			if ( seqpos == CLAIMED_RESID ) {
				TS_ASSERT_DELTA( final_pose.phi( seqpos ), NEW_PHI, 0.000001 );
			} else {
				TS_ASSERT_DELTA( final_pose.phi( seqpos ), init_phis[ seqpos ], 0.000001 );
			}
		}

		TR <<  "End: test_dual_environment"  << std::endl;
	}

	void test_single_phi_moves(){
		TR <<  "Beginning: test_single_phi_moves"  << std::endl;
		using namespace protocols::environment;
		using namespace core::environment;

		//Tests that
		TorsionMoverOP allowed_mover( new TorsionMover( true, true ) );
		TorsionMoverOP duplicate_claim_mover( new TorsionMover( true, true ) );
		TorsionMoverOP no_claim_mover( new TorsionMover( false, true ) );
		TorsionMoverOP unreg_mover( new TorsionMover( true, true ) );

		EnvironmentOP env_op( new Environment( "torsion" ) );
		Environment & env = *env_op;

		env.register_mover( allowed_mover );
		env.register_mover( duplicate_claim_mover );
		env.register_mover( no_claim_mover );
		//don't register unreg_mover

		core::pose::Pose final_pose;

		{
			core::pose::Pose protected_pose = env.start( pose );
			// Verify conformation got copied into protected_pose pose.
			TS_ASSERT_EQUALS( protected_pose.size(), pose.size() );

			// Verify no_claim_mover can't change anything -- it shouldn't have a passport for this environment (NullPointer excn)
			TS_ASSERT_THROWS( no_claim_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
			TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

			// Verify no_lock_mover can't change anything -- protected_pose shouldn't have a passport on its unlock stack
			TS_ASSERT_THROWS( allowed_mover->missing_unlock_apply( protected_pose ), EXCN_Env_Security_Exception );
			TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

			// Verify that unregistered mover lacks a passport for protected_pose conformation
			TS_ASSERT_THROWS( unreg_mover->apply( protected_pose ), utility::excn::NullPointerError );
			TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

			// Verify that allowed_mover can change it's claimed angle
			TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
			TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), NEW_PHI, 1e-12 );

			// Verify that allowed_mover can do it again.
			TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
			TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), NEW_PHI, 1e-10 );

			//Verify that duplicate mover *also* can make the change without throwing an error
			TS_ASSERT_THROWS_NOTHING( duplicate_claim_mover->apply( protected_pose ) );
			//Verify that allowed can't change 9 while claiming CLAIMED_RESID
			TS_ASSERT_THROWS( allowed_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
			TS_ASSERT_THROWS( duplicate_claim_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
			TS_ASSERT_DELTA( protected_pose.phi( UNCLAIMED_RESID ), init_phis[ UNCLAIMED_RESID ], 1e-12 );

			//Verify angles 1-9 are untouched in protected_pose
			for ( core::Size i = 1; i <= pose.size()-1; ++i ) {
				TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 1e-12 );
			}

			TS_ASSERT_THROWS_NOTHING( final_pose = env.end( protected_pose ) );
		}

		// Verify angles 1-9 are untouched in pose and final_pose;
		// Phi 1 is not well-defined, so skip that one.
		for ( core::Size i = 2; i <= pose.size(); ++i ) {
			if ( i != CLAIMED_RESID ) {
				TS_ASSERT_DELTA( pose.phi( i ), final_pose.phi( i ), 1e-12 );
				TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 1e-12 );
			}
		}

		//Finish verification inital pose is unaffected
		TS_ASSERT_DELTA( pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

		//Verify angle changes propagate out of environment
		TS_ASSERT_DELTA( final_pose.phi( CLAIMED_RESID ), NEW_PHI, 1e-12 );

		//Verify other angle changes don't back-propagate to original pose
		TS_ASSERT_DIFFERS( pose.phi( CLAIMED_RESID ), final_pose.phi( CLAIMED_RESID ) );

		TR <<  "End: test_single_phi_moves"  << std::endl;
	}

	void test_torsion_must_can_coexist(){
		TR <<  "Beginning: test_torsion_must_can_coexist"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		// Veirfy coexeistance of must and can control
		TorsionMoverOP must_mover( new TorsionMover( true, true, claims::MUST_CONTROL ) );
		TorsionMoverOP can_mover( new TorsionMover( true, true, claims::CAN_CONTROL ) );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;

		env.register_mover( must_mover );
		env.register_mover( can_mover );

		core::pose::Pose protected_pose;
		TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
		TS_ASSERT_THROWS_NOTHING( must_mover->apply( protected_pose ) );
		TS_ASSERT_THROWS_NOTHING( can_mover->apply( protected_pose ) );
		TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );

		TR <<  "End: test_torsion_must_can_coexist"  << std::endl;
	}

	void test_torsion_can_exclusive_compatibility(){
		TR <<  "Beginning: torsion_can_exclusive_compatibility"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		TorsionMoverOP exclusive_mover( new TorsionMover( true, true, claims::EXCLUSIVE ) );
		TorsionMoverOP can_mover( new TorsionMover( true, true, claims::CAN_CONTROL ) );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;

		env.register_mover( exclusive_mover );
		env.register_mover( can_mover );

		core::pose::Pose protected_pose;
		TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
		TS_ASSERT_THROWS( can_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
		TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ] );

		TS_ASSERT_THROWS_NOTHING( exclusive_mover->apply( protected_pose ) );
		TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );

		TR <<  "End: torsion_can_exclusive_compatibility"  << std::endl;
	}

	void test_torsion_must_exclusive_incompatibility(){
		TR <<  "Beginning: test_torsion_must_exclusive_incompatibility"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		TorsionMoverOP exclusive_mover( new TorsionMover( true, true, claims::EXCLUSIVE ) );
		TorsionMoverOP must_mover( new TorsionMover( true, true, claims::MUST_CONTROL ) );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;

		env.register_mover( exclusive_mover );
		env.register_mover( must_mover );

		core::pose::Pose protected_pose;
		TS_ASSERT_THROWS( protected_pose = env.start( pose ), utility::excn::BadInput );

		TR <<  "End: test_torsion_must_exclusive_incompatibility"  << std::endl;
	}

	void test_torsion_init(){
		TR <<  "Beginning: test_torsion_init"  << std::endl;

		using namespace protocols::environment;
		using namespace core::environment;

		TorsionMoverOP must_init_mover( new TorsionMover( true, true,
			claims::DOES_NOT_CONTROL,
			claims::EXCLUSIVE ) );
		TorsionMoverOP must_init_mover2( new TorsionMover( true, true,
			claims::DOES_NOT_CONTROL,
			claims::EXCLUSIVE ) );
		TorsionMoverOP active_can_init_mover( new TorsionMover( true, true,
			claims::DOES_NOT_CONTROL,
			claims::CAN_CONTROL ) );
		TorsionMoverOP inactive_can_init_mover( new TorsionMover( true, false,
			claims::DOES_NOT_CONTROL,
			claims::CAN_CONTROL ) );

		core::pose::Pose protected_pose;

		{ // 2xMUST_INIT are not compatible
			EnvironmentOP env_op( new Environment( "env" ) );
			Environment & env = *env_op;

			env.register_mover( must_init_mover );
			env.register_mover( must_init_mover2 );
			TS_ASSERT_THROWS( protected_pose = env.start( pose ), utility::excn::BadInput );
		}

		{ // MUST_INIT takes precedence over CAN_INIT and does not complain.
			EnvironmentOP env_op( new Environment( "env" ) );
			Environment & env = *env_op;

			env.register_mover( must_init_mover );
			env.register_mover( inactive_can_init_mover );
			TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
			TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );
		}

		TR <<  "End: test_torsion_init"  << std::endl;
	}

};
