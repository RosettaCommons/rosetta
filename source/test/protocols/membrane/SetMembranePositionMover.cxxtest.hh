// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/SetMembranePositionMover.cxxtest.hh
///
/// @brief		Unit Test for Membrane Position Move (Rotation & Translation, Uniform)
/// @details	The memrbane fold tree assumes that the membrane residue is attached by jump to
///				the rest of the pose. This assumption asserts that all positions explicitly describing
///				the membrane will be downstream of this jump. This code should apply the correct translation
///				and rotation to describe a new position of the membrane.
///				Last Updated: 7/8/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

// C++ Headers
#include <cstdlib>

using namespace core;

/// @brief Test Suite for Membrane Embedding factory
class SetMembranePositionTest : public CxxTest::TestSuite {

public:

    /// @brief Setup
    void setUp()
    {

		using namespace basic::options;
		using namespace core::conformation::membrane;
		using namespace protocols::membrane;

		// Initialize Rosetta
        protocols_init();

		// Load Pose from pdb
		std::string input_pose = "protocols/membrane/1C3W_TR_A.pdb";
		pose_ = core::import_pose::pose_from_pdb( input_pose );

		// Load span object from spanfile
		std::string spanfile = "protocols/membrane/1C3W_A.span";

		// Define a starting membrane position
		Vector center( 0, 0, 0 );
		Vector normal( 0, 0, 15 );

		// Add Membrane to pose!
		AddMembraneMoverOP add_memb( new AddMembraneMover( center, normal, spanfile, 1 ) );
		add_memb->apply( *pose_ );

    }

    /// @brief teardown
    void tearDown()
    {}

    /// @brief Testing Consecutive Transformations of the Membrane Center
    void test_membrane_position_translation() {

		TS_TRACE( "Testing Deterministic Translation of the Membrane Position" );

		using namespace protocols::membrane;

		// Translate back along all 3 axes
		Vector xyz_trans( -5, -11, -15 );

		// Apply translation move
		SetMembraneCenterMoverOP xyz_move( new SetMembraneCenterMover( xyz_trans ) );
		xyz_move->apply( *pose_ );

		// Check result
		TS_ASSERT( xyz_trans == pose_->conformation().membrane_info()->membrane_center() );

    }

    /// @brief Testing Consecutive Rotations of the membrane Normal (about fixed center)
    void test_membrane_position_rotation() {

		TS_TRACE( "Testing Deterministic Rotation of the membrane position" );

		using namespace protocols::membrane;

		// Grab current normal
		Vector current_normal( pose_->conformation().membrane_info()->membrane_normal() );

		// Simple rotation 1
		Vector rot_1( 0, 15, 0 );
		SetMembraneNormalMoverOP first_move( new SetMembraneNormalMover( rot_1 ) );
		first_move->apply( *pose_ );
		TS_ASSERT_DELTA( angle_of( current_normal, pose_->conformation().membrane_info()->membrane_normal() ), 1.57, 0.001);
		TS_ASSERT( position_equal_within_delta( rot_1, pose_->conformation().membrane_info()->membrane_normal(), 0.0001 ) );

	}

	/// @brief Testing uniform rotation & translation
	void test_uniform_rotation_translation() {

		TS_TRACE( "Testing rotation & translation move" );

		using namespace protocols::membrane;

		// Pick a new center/normal position
		Vector new_center( 0, 5, 10 );
		Vector new_normal( 0, 15, 0 );

		// Apply Rotation and translation move
		SetMembranePositionMoverOP rt( new SetMembranePositionMover( new_center, new_normal ) );
		rt->apply( *pose_ );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( new_center, pose_->conformation().membrane_info()->membrane_center(), 0.0001 ) );
		TS_ASSERT( position_equal_within_delta( new_normal, pose_->conformation().membrane_info()->membrane_normal(), 0.0001 ) );
	}

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}


private: // data

    // Resulting Membrane Protein
    core::pose::PoseOP pose_;

}; // class SetMembranePositionTest

