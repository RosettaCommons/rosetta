// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/TransformIntoMembraneMover.cxxtest.hh
/// @brief		Unit Test for transforming pose into fixed membrane
/// @details	Testing for correct translation and rotation of the
///				pose into the membrane coordinate frame
/// @author		JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/geometry/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/conformation/membrane/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

// C++ Headers
#include <cstdlib>

using namespace core;

/// @brief Test Suite for transformin a pose into membrane coordinates
class TransformIntoMembraneMoverTest : public CxxTest::TestSuite {

public:

    /// @brief Setup
    void setUp()
    {
		using namespace basic::options;
		using namespace core::conformation::membrane;
		using namespace protocols::membrane::geometry;
		using namespace protocols::membrane;

		// Initialize Rosetta
        protocols_init();

		// load pose
		std::string pdbfile = "protocols/membrane/1AFO_AB_before_out.pdb";
		pose_ = core::import_pose::pose_from_pdb( pdbfile );

		// Load span object from spanfile
		spanfile_ = "protocols/membrane/1AFO_AB.span";

		// Add Membrane to pose
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile_, 0 ) );
		add_memb->apply( *pose_ );

		// reorder foldtree
		reorder_membrane_foldtree( *pose_ );

    }

    /// @brief teardown
    void tearDown()
    {}

	/// @brief test transform into default membrane
	void test_transform_into_default_membrane() {

		TS_TRACE( "TESTING TRANSFORM INTO DEFAULT MEMBRANE" );

		using namespace protocols::membrane;

		// set new membrane to transform pose into
		Vector new_center ( mem_center );
		Vector new_normal ( mem_normal );

		// membrane default that are not moved, output in PDB
		Vector m_center ( mem_center );
		Vector m_normal ( mem_normal );
		Vector m_thickness( mem_thickness, 0, 0 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( new_center, new_normal, spanfile_ ) );
		transform->apply( *pose_ );

		// check positions of CA atoms of first and last residue after rotation
		Vector res1_after ( -11.5162, 14.2177, -3.9571 );
		Vector res40_after (  -14.1853, -4.9844, 24.6974 );
		Vector res41_after ( 9.1834, 1.7434, -5.6530 );
		Vector res80_after ( 11.2935, -6.9151, 23.5296 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

			// check positions of center and normal
		TS_ASSERT( position_equal_within_delta( m_thickness, pose_->residue(81).atom(1).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->residue(81).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->residue(81).atom(3).xyz(), 0.001 ) );

	}

	/// @brief test transform into user-defined membrane
	void test_transform_into_userdefined_membrane() {

		TS_TRACE( "TESTING TRANSFORM INTO USER-DEFINED MEMBRANE" );

		using namespace protocols::membrane;

		// set new membrane to transform pose into
		Vector new_center (20, 20, 20);
		Vector new_normal (15, 0, 0);

		// membrane default that are not moved, output in PDB
		Vector m_center ( mem_center );
		Vector m_normal ( mem_normal );
		Vector m_thickness( mem_thickness, 0, 0 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( new_center, new_normal, spanfile_ ) );
		transform->apply( *pose_ );

		// check positions of CA atoms of first and last residue after rotation
		Vector res1_after ( 16.0428, 37.8056, 24.2101 );
		Vector res40_after ( 44.6974, 21.6612, 34.9435 );
		Vector res41_after ( 14.3469, 17.5878, 10.9690 );
		Vector res80_after ( 43.5296, 8.8708, 12.8232 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT( position_equal_within_delta( m_thickness, pose_->residue(81).atom(1).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->residue(81).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->residue(81).atom(3).xyz(), 0.001 ) );
	}

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}


private: // data

    // store data
    core::pose::PoseOP pose_;
	std::string spanfile_;

}; // class TransformIntoMembraneMoverTest

