// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/TransformIntoMembraneMover.cxxtest.hh
/// @brief  Unit Test for transforming pose into fixed membrane
/// @details Testing for correct translation and rotation of the
///    pose into the membrane coordinate frame
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project Headers
#include <protocols/membrane/util.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>


// Package Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

using namespace core;
using namespace core::conformation::membrane;
using namespace protocols::membrane::geometry;
using namespace protocols::membrane;

/// @brief Test Suite for transformin a pose into membrane coordinates
class TransformIntoMembraneMoverTest : public CxxTest::TestSuite {

public:

	/// @brief Setup
	void setUp()
	{
		using namespace core::conformation::membrane;
		using namespace protocols::membrane::geometry;
		using namespace protocols::membrane;

		// Initialize Rosetta
		protocols_init();

		// Test coordinate transformation applied to individual residues is correct
		// Case: Glycophorin A
		std::string pdbfile = "protocols/membrane/1AFO_AB_before_out.pdb";
		pose_ = core::import_pose::pose_from_pdb( pdbfile );
		Vector m_center( 0, 0, 0 );
		Vector m_normal( 0, 0, 1 );
		AddMembraneMoverOP add_memb( new AddMembraneMover( m_center, m_normal, "protocols/membrane/1AFO_AB.span", 0 ) );
		add_memb->apply( *pose_ );

		// Test angle between single helix and normal azis is 0 after transformation
		// Case: TM domain of the M2 proton channel (single helix)

		m2_pose_ = PoseOP( new Pose() );
		core::import_pose::pose_from_pdb( *m2_pose_, "protocols/membrane/1mp6.pdb" );
		AddMembraneMoverOP add_memb1 = AddMembraneMoverOP( new AddMembraneMover( "protocols/membrane/1mp6.span" ) );
		add_memb1->apply( *m2_pose_ );

	}

	/// @brief teardown
	void tearDown()
	{}

	////////////////////////////////////////////////////////////////////////////

	/// @brief test transform into default membrane
	void test_transform_into_default_membrane() {

		TS_TRACE( "=======Testing transform into membrane mover on default membrane coordinates" );

		using namespace protocols::membrane;

		// membrane default that are not moved, output in PDB
		Vector m_center ( 0, 0, 0 );
		Vector m_normal ( 0, 0, 1 );
		Real m_thickness( 15 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
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
		TS_ASSERT_DELTA( m_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->conformation().membrane_info()->membrane_center(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->conformation().membrane_info()->membrane_normal(), 0.001 ) );

	}

	////////////////////////////////////////////////////////////////////////////

	/// TODO: test constructor from jump

	////////////////////////////////////////////////////////////////////////////

	// test constructor from embedding
	void test_constructor_from_embedding() {

		TS_TRACE( "=======Testing transform into membrane mover from defined embedding" );

		using namespace protocols::membrane;

		// set the current embedding of the protein
		Vector cur_emb_center (20.9885, -0.04375, 56.6125);
		Vector cur_emb_normal (-8.6212, 0.83567, 17.05025);
		EmbeddingDefOP current_emb( new EmbeddingDef( cur_emb_center, cur_emb_normal ) );

		// membrane default that are not moved, output in PDB
		Vector m_center ( 0, 0, 0 );
		Vector m_normal ( 0, 0, 1 );
		Real m_thickness( 15 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( current_emb ) );
		transform->apply( *pose_ );

		// check positions of CA atoms of first and last residue after rotation
		Vector res1_after ( -8.317, 10.922, -16.568 );
		Vector res40_after ( -22.831, -3.756, 11.196 );
		Vector res41_after ( 11.119, -0.315, -7.461 );
		Vector res80_after ( 0.795, -4.211, 20.987 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT_DELTA( m_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->conformation().membrane_info()->membrane_center(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->conformation().membrane_info()->membrane_normal(), 0.001 ) );
	}

	////////////////////////////////////////////////////////////////////////////

	/// @brief test transform into user-defined membrane
	void test_transform_into_userdefined_membrane() {

		TS_TRACE( "=======Testing transform into membrane mover on user provided membrane coordinates" );

		using namespace protocols::membrane;

		// set new membrane to transform pose into
		Vector m_center (20, 20, 20);
		Vector m_normal (1, 0, 0);
		Real m_thickness( 15 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( m_center, m_normal ) );
		transform->apply( *pose_ );

		// check positions of CA atoms of first and last residue after rotation
		Vector res1_after ( 16.043, 37.806, 24.210 );
		Vector res40_after ( 44.697, 21.661, 34.944 );
		Vector res41_after ( 14.347, 17.588, 10.969 );
		Vector res80_after ( 43.530, 8.871, 12.823 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT_DELTA( m_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->conformation().membrane_info()->membrane_center(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->conformation().membrane_info()->membrane_normal(), 0.001 ) );
	}

	////////////////////////////////////////////////////////////////////////////

	// test constructor form embedding and user-defined membrane
	void test_constructor_from_embedding_and_membrane() {

		TS_TRACE( "=======Testing transform into membrane mover from embedding and defined membrane" );

		using namespace protocols::membrane;

		// set the current embedding of the protein
		Vector cur_emb_center (20.9885, -0.04375, 56.6125);
		Vector cur_emb_normal (-8.6212, 0.83567, 17.05025);
		EmbeddingDefOP current_emb( new EmbeddingDef( cur_emb_center, cur_emb_normal ) );

		// set new membrane to transform pose into
		Vector m_center (20, 20, 20);
		Vector m_normal (1, 0, 0);
		Real m_thickness( 15 );

		// Apply Rotation and translation move
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( current_emb, m_center, m_normal ) );
		transform->apply( *pose_ );

		// check positions of CA atoms of first and last residue after rotation
		Vector res1_after ( 3.432, 31.406, 27.639 );
		Vector res40_after ( 31.196, 17.635, 43.017 );
		Vector res41_after ( 12.539, 19.011, 8.921 );
		Vector res80_after ( 40.987, 15.749, 19.462 );

		// Check the structure was moved to the correct position
		TS_ASSERT( position_equal_within_delta( res1_after, pose_->residue(1).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res40_after, pose_->residue(40).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res41_after, pose_->residue(41).atom(2).xyz(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( res80_after, pose_->residue(80).atom(2).xyz(), 0.001 ) );

		// check positions of center and normal
		TS_ASSERT_DELTA( m_thickness, pose_->conformation().membrane_info()->membrane_thickness(), 0.001 );
		TS_ASSERT( position_equal_within_delta( m_center, pose_->conformation().membrane_info()->membrane_center(), 0.001 ) );
		TS_ASSERT( position_equal_within_delta( m_normal, pose_->conformation().membrane_info()->membrane_normal(), 0.001 ) );
	}

	////////////////////////////////////////////////////////////////////////////

	/// @brief Check that when a protein is transformed, the angle between the protein "axis" and
	/// membrane normal axis is 0 (the "definition" of correct transformation)
	void test_helix_to_normal_angle() {

		TS_TRACE( "=======Testing helix to normal angle" );

		// Apply transformation to the protein
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
		transform->apply( *m2_pose_ );

		// Get membrane normal from membrane info
		Vector normal( m2_pose_->conformation().membrane_info()->membrane_normal() );

		// Calculate the "vector" representing the helix as the difference
		// between start & end TM helices.
		SpanOP helix_span( m2_pose_->conformation().membrane_info()->spanning_topology()->span( 1 ) );
		EmbeddingDefOP span_embed = EmbeddingDefOP( new EmbeddingDef( *m2_pose_, helix_span->start(), helix_span->end() ) );
		Vector helix_axis( span_embed->normal() );
		normal.normalize();
		helix_axis.normalize();

		// Check that both vectors have a positive orientation
		if ( normal.length() > 0 ) normal.negate();
		if ( helix_axis.length() > 0 ) helix_axis.negate();

		// Check tilt angle is equal to zero
		Real tilt_angle( angle_of( normal, helix_axis ) );
		TS_ASSERT_DELTA( tilt_angle, 0.0, 0.005 );

	}

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}


private: // data

	// Per-residue testing
	core::pose::PoseOP pose_;

	// Helix angles testing
	core::pose::PoseOP m2_pose_;

}; // class TransformIntoMembraneMoverTest

