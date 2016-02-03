// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/TranslationRotationMover.cxxtest.hh
/// @brief  Unit Test for translating and rotating pose into membrane
/// @details Testing for correct translation and rotation of the
///    pose into the membrane coordinate frame
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// Utility
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>

using namespace core;
using namespace utility;

/// @brief Test Suite for transformin a pose into membrane coordinates
class TranslationRotationMoverTest : public CxxTest::TestSuite {

public:

	/// @brief Setup
	void setUp()
	{
		using namespace basic::options;
		using namespace core::conformation::membrane;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;

		// Initialize Rosetta
		protocols_init();

		// load pose
		pose_ = core::import_pose::pose_from_file( "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);
		std::string spanfile = "protocols/membrane/1AFO_AB.span";

		// Add Membrane to pose
		AddMembraneMoverOP add_memb( new AddMembraneMover( spanfile, 1 ) );
		add_memb->apply( *pose_ );

		// reorder the foldtree to have membrane residue at the root
		reorder_membrane_foldtree( *pose_ );

		// get membrane jump
		memjump_ = pose_->conformation().membrane_info()->membrane_jump();
		TS_TRACE( "mem residue: " + to_string(pose_->conformation().membrane_info()->membrane_rsd_num()) );
		TS_TRACE( "mem jump: " + to_string(pose_->conformation().membrane_info()->membrane_jump()));
		pose_->fold_tree().show( std::cout );

	}

	/// @brief teardown
	void tearDown()
	{}

	/// @brief test translation of membrane pose
	void test_translation() {

		TS_TRACE( "TESTING TRANSLATION MOVE" );

		using namespace core::pose;
		using namespace protocols::membrane;

		// define vectors
		Vector trans(10, 10, 10);

		// Apply Rotation and translation move
		TranslationMoverOP translate( new TranslationMover( trans, memjump_ ) );
		TS_TRACE( "mem residue: " + to_string(pose_->conformation().membrane_info()->membrane_rsd_num()) );
		TS_TRACE( "mem jump: " + to_string(pose_->conformation().membrane_info()->membrane_jump()));
		pose_->fold_tree().show( std::cout );
		translate->apply( *pose_ );

		// check positions of CA atoms of first and last residue after translation
		Vector res1_after (-7.655, 5.184, 5.480);
		Vector res40_after (6.425, 0.503, 36.734);
		Vector res41_after (12.974, 17.504, 2.374);
		Vector res80_after (19.846, 21.712, 31.804);

		// check positions of new center and normal
		Vector new_center ( 0, 0, 0 );
		Vector new_normal ( 0, 0, 1 );
		Real new_thickness ( 15 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT_DELTA( new_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( new_center, pose_->conformation().membrane_info()->membrane_center( pose_->conformation() ), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( new_normal, pose_->conformation().membrane_info()->membrane_normal( pose_->conformation() ), 0.001 ) );
	}

	/// @brief test rotation of membrane pose
	void test_rotation() {

		TS_TRACE( "TESTING ROTATION MOVE" );

		using namespace core::pose;
		using namespace protocols::membrane;

		// define vectors
		Vector old_normal(0, 0, 1);
		Vector new_normal(1, 0, 0);
		Vector rot_center(0, 0, 0);

		// Apply Rotation and translation move
		RotationMoverOP rotate( new RotationMover( old_normal, new_normal, rot_center, memjump_ ) );
		TS_TRACE( "mem residue: " + to_string(pose_->conformation().membrane_info()->membrane_rsd_num()) );
		TS_TRACE( "mem jump: " + to_string(pose_->conformation().membrane_info()->membrane_jump()));
		pose_->fold_tree().show( std::cout );
		rotate->apply( *pose_ );

		// check positions of CA atoms of first and last residue after translation
		Vector res1_after (-4.520, -4.816, 17.655);
		Vector res40_after (26.734, -9.497, 3.575);
		Vector res41_after (-7.626, 7.504, -2.974);
		Vector res80_after (21.804, 11.712, -9.846);

		// check that membrane didn't move
		Vector m_center ( 0, 0, 0 );
		Vector m_normal ( 0, 0, 1 );
		Real m_thickness ( 15 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT_DELTA( m_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->conformation().membrane_info()->membrane_center( pose_->conformation() ), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->conformation().membrane_info()->membrane_normal( pose_->conformation() ), 0.001 ) );
	}

	/// @brief test rotation and translation of membrane pose
	void test_rotation_translation() {

		TS_TRACE( "TESTING TRANSLATION AND ROTATION MOVES" );

		using namespace core::pose;
		using namespace protocols::membrane;

		// define vectors
		Vector old_center(0, 0, 0);
		Vector old_normal(0, 0, 1);
		Vector new_center(20, 20, 20);
		Vector new_normal(1, 0, 0);

		// Apply Rotation and translation move
		TranslationRotationMoverOP rt( new TranslationRotationMover( old_center, old_normal, new_center, new_normal, memjump_ ) );
		TS_TRACE( "mem residue: " + to_string(pose_->conformation().membrane_info()->membrane_rsd_num()) );
		TS_TRACE( "mem jump: " + to_string(pose_->conformation().membrane_info()->membrane_jump()));
		pose_->fold_tree().show( std::cout );
		rt->apply( *pose_ );

		// check positions of CA atoms of first and last residue after translation
		Vector res1_after (15.480, 15.184, 37.655);
		Vector res40_after (46.734, 10.503, 23.575);
		Vector res41_after (12.374, 27.504, 17.026);
		Vector res80_after (41.804, 31.712, 10.154);

		// check that membrane didn't move
		Vector m_center ( 0, 0, 0 );
		Vector m_normal ( 0, 0, 1 );
		Vector m_thickness ( 15, 0, 0);

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
	Size memjump_;

}; // class TranslationRotationMoverTest

