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

#include <basic/Tracer.hh>

#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/rigid/UniformRigidBodyCM.hh>

#include <test/core/init_util.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

//C++ headers
#include <iostream>

using namespace core;

class UniformRigidBody : public CxxTest::TestSuite {
public:

	// Shared data elements go here.
	core::pose::Pose pose;

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();

		core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			pose.set_phi( i, -65 );
			pose.set_psi( i, -41 );
		}

		pose.append_pose_by_jump( *pose.clone(), 5 );
		core::kinematics::Jump j = pose.jump( 1 );
		j.gaussian_move(1, 50, 0);
		pose.set_jump(1, j);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_uniform_rb() {

		using namespace core::environment;
		using namespace protocols;
		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace numeric;
		using namespace core::conformation;
		//using namespace abinitio::abscript;
		using namespace protocols::rigid;

		TS_ASSERT_DIFFERS( pose.total_residue(), 0 );

		UniformRigidBodyCMOP rigpert( new UniformRigidBodyCM( "perturb", LocalPosition( "BASE", 1 ), LocalPosition( "BASE", pose.total_residue() ) ) );
		UniformRigidBodyCMOP rigpert_dup( new UniformRigidBodyCM( "perturb2", LocalPosition( "BASE", 1 ), LocalPosition( "BASE", pose.total_residue() ) ) );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;
		env.register_mover( rigpert );
		env.register_mover( rigpert_dup );

		core::pose::Pose ppose;

		//Test correct topology construction
		TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
		TS_ASSERT_EQUALS( ppose.num_jump(), 1 );
		TS_ASSERT( ppose.fold_tree().is_jump_point( 1 ) );
		TS_ASSERT( ppose.fold_tree().is_jump_point( ppose.total_residue() ) );

		//Verify initialization perturbation
		core::kinematics::Jump j;
		{
			core::pose::Pose tmp( pose );
			tmp.fold_tree( ppose.fold_tree() );
			j = tmp.jump( 1 );
		}

		TS_ASSERT_DIFFERS( j.get_translation(), ppose.jump( 1 ).get_translation() );
		TS_ASSERT_DIFFERS( j.get_rotation(), ppose.jump( 1 ).get_rotation() );
		j = ppose.jump( 1 );

		//Verify changes after perturbation.
		TS_ASSERT_THROWS_NOTHING( rigpert->apply( ppose ) );
		TS_ASSERT_DIFFERS( j.get_translation(), ppose.jump( 1 ).get_translation() );
		TS_ASSERT_DIFFERS( j.get_rotation(), ppose.jump( 1 ).get_rotation() );

		// Verify that duplicate movers can share.
		TS_ASSERT_THROWS_NOTHING( rigpert_dup->apply( pose ) );

		//Test environment end
		core::pose::Pose final_pose;
		TS_ASSERT_THROWS_NOTHING( final_pose = env.end( ppose ) );
	}

};
