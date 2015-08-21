// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>
#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/Environment.hh>

#include <test/protocols/environment/TestClientMover.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

//C++ headers
#include <iostream>
#include <boost/bind/bind.hpp>

// --------------- Test Class --------------- //

class ClaimTests : public CxxTest::TestSuite {
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
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			init_phis.push_back( pose.phi( i ) );
		}

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_tclaim_constructor_single_pos() {

		using namespace protocols::environment;
		using namespace protocols::environment::claims;

		core::Size const SEQPOS = 5;

		protocols::environment::EnvironmentOP env_op( new protocols::environment::Environment( "env" ) );
		protocols::environment::Environment & env = *env_op;

		TorsionClaimOP claim( new TorsionClaim( NULL, core::environment::LocalPosition( "BASE", SEQPOS ) ) );
		claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
		claim->claim_sidechain( false );
		claim->claim_backbone( true );
		TesterOP claim_test( new Tester );
		claim_test->init( claim );

		env.register_mover( claim_test );

		core::pose::Pose prot_pose;
		TS_ASSERT_THROWS_NOTHING( prot_pose = env.start( pose ) );

		void (core::pose::Pose::*set_phi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_phi;
		void (core::pose::Pose::*set_psi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_psi;
		void (core::pose::Pose::*set_omega) ( core::Size const, core::Real const ) = &core::pose::Pose::set_omega;

		for ( core::Size i = 1; i <= prot_pose.total_residue(); ++i ) {
			if ( i != SEQPOS ) {
				std::cout << i << std::endl;
				if ( i != 1 ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
				if ( i != prot_pose.total_residue() ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
			}
		}

		// Exception safety of the above calls is tested in other tests.

		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, SEQPOS, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, SEQPOS, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, SEQPOS, 0.0 ) ) );
	}

	void test_tclaim_constructor_multi_pos() {

		using namespace core::environment;
		using namespace protocols::environment;
		using namespace protocols::environment::claims;

		core::Size const SEQPOS1 = 5;
		core::Size const SEQPOS2 = 7;

		protocols::environment::EnvironmentOP env_op( new protocols::environment::Environment( "env" ) );
		protocols::environment::Environment & env = *env_op;

		LocalPositions envpos;
		envpos.push_back( core::environment::LocalPositionOP( new LocalPosition( "BASE", SEQPOS1 ) ) );
		envpos.push_back( core::environment::LocalPositionOP( new LocalPosition( "BASE", SEQPOS2 ) ) );

		TorsionClaimOP claim( new TorsionClaim( NULL, envpos ) );
		claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
		claim->claim_sidechain( false );
		claim->claim_backbone( true );
		TesterOP claim_test( new Tester() );
		claim_test->init( claim );

		env.register_mover( claim_test );

		core::pose::Pose prot_pose;
		TS_ASSERT_THROWS_NOTHING( prot_pose = env.start( pose ) );

		void (core::pose::Pose::*set_phi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_phi;
		void (core::pose::Pose::*set_psi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_psi;
		void (core::pose::Pose::*set_omega) ( core::Size const, core::Real const ) = &core::pose::Pose::set_omega;

		for ( core::Size i = 1; i <= prot_pose.total_residue(); ++i ) {
			if ( i != SEQPOS1 && i != SEQPOS2 ) {
				std::cout << i << std::endl;
				if ( i != 1 ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
				if ( i != prot_pose.total_residue() ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
			}
		}

		// Exception safety of the above calls is tested in other tests.

		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, SEQPOS1, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, SEQPOS1, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, SEQPOS1, 0.0 ) ) );

		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, SEQPOS2, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, SEQPOS2, 0.0 ) ) );
		TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, SEQPOS2, 0.0 ) ) );
	}

	void test_tclaim_constructor_list_pos( core::pose::Pose const& pose ) {

		using namespace core::environment;
		using namespace protocols::environment;
		using namespace protocols::environment::claims;

		core::Size const SEQPOS_START = 4;
		core::Size const SEQPOS_END = 7;

		protocols::environment::EnvironmentOP env_op( new protocols::environment::Environment( "env" ) );
		protocols::environment::Environment & env = *env_op;

		TorsionClaimOP claim( new TorsionClaim( NULL, "BASE", std::make_pair( SEQPOS_START, SEQPOS_END ) ) );
		claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
		claim->claim_sidechain( false );
		claim->claim_backbone( true );
		TesterOP claim_test( new Tester() );
		claim_test->init( claim );

		env.register_mover( claim_test );

		core::pose::Pose prot_pose;
		TS_ASSERT_THROWS_NOTHING( prot_pose = env.start( pose ) );

		void (core::pose::Pose::*set_phi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_phi;
		void (core::pose::Pose::*set_psi)   ( core::Size const, core::Real const ) = &core::pose::Pose::set_psi;
		void (core::pose::Pose::*set_omega) ( core::Size const, core::Real const ) = &core::pose::Pose::set_omega;

		for ( core::Size i = 1; i <= prot_pose.total_residue(); ++i ) {
			if ( i < SEQPOS_START || i > SEQPOS_END ) {
				std::cout << i << std::endl;
				if ( i != 1 ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
				if ( i != prot_pose.total_residue() ) {
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
					TS_ASSERT_THROWS( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, i, 0.0 ) ) , EXCN_Env_Security_Exception );
				}
			}
		}

		// Exception safety of the above calls is tested in other tests.

		for ( core::Size i = SEQPOS_START; i <= SEQPOS_END; ++i ) {
			TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_phi, &prot_pose, i, 0.0 ) ) );
			TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_psi, &prot_pose, i, 0.0 ) ) );
			TS_ASSERT_THROWS_NOTHING( claim_test->apply( prot_pose, boost::bind( set_omega, &prot_pose, i, 0.0 ) ) );
		}
	}
};
