// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/protocols/carbohydrates/RingPlaneFlipMover.cxxtest.hh
/// @brief   Test suite for the RingPlaneFlipMover
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <protocols/carbohydrates/RingPlaneFlipMover.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Numeric header
#include <numeric/angle.functions.hh>

// Basic header
#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.carbohydrates.RingPlaneFlipMover.cxxtest" );


class RingPlaneFlipMoverTests : public CxxTest::TestSuite {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core;
		using namespace pose;

		core_init_with_additional_options( "-include_sugars" );

		pose_ = PoseOP( new Pose );
		make_pose_from_saccharide_sequence( *pose_,
			"b-D-Glcp-(1->4)-[b-D-Glcp-(1->6)]-b-D-Glcp-(1->4)-b-D-Manp-(1->4)-b-D-Glcp-(1->4)-"
			"b-D-Galp-(1->4)-a-D-Glcp-(1->3)-b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp" );
		// TODO: Replace residue 7 with a 5-linked, all-equatorial ketose.
		extend_pose( *pose_ );
	}

	// Destruction
	void tearDown()
	{}

	// Tests ///////////////////////////////////////////////////////////////////
	// Confirm that the proper number of conformers are loaded into 5- and 6-membered ring sets.
	void test_RingPlaneFlipMover()
	{
		using namespace std;
		using namespace numeric;
		using namespace core::kinematics;
		using namespace core::select::residue_selector;

		protocols::carbohydrates::RingPlaneFlipMover mover;
		ResidueIndexSelectorOP selector( new ResidueIndexSelector );
		MoveMapOP mm( new MoveMap );
		mm->set_bb( true );
		mm->set_chi( false );
		mm->set_branches( true );
		mover.movemap( mm );

		TR << "Testing that RingPlaneFlipMover only flips the rings that it should." << endl;

		TR << " Residue 1: flip (lower terminus; only psi2 should move)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "1" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 2: flip (4-eq, 4-linked, & beta)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "2" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 3: no flip (3-linked)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "3" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 4: no flip (alpha)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "4" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 5: no flip (Gal is 4-ax)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "5" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 6: flip (5-eq, 5-linked, & beta ketose)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "6" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 7: flip (4-eq, (2-ax,) 4-linked, & beta ketose)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "7" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 8: no flip (branch point)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "8" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 9: flip  (upper terminus; only phi9 should move)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "9" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ),   0, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ), 180, 1 );

		TR << " Residue 10: flip (upper terminus; only phi10 should move)" << endl;
		extend_pose( *pose_ );
		selector->set_index( "10" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  2 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  3 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  4 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  5 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  6 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  7 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  8 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi(  9 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->psi( 10 ) ) ), 180, 1 );
		TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( int( pose_->phi( 10 ) ) ),   0, 1 );
	}

private:
	// Private function ////////////////////////////////////////////////////////
	void extend_pose( core::pose::Pose & pose )
	{
		for ( core::uint i( 1 ); i <= pose.size(); ++i ) {
			pose.set_phi( i, 180.0 );
			pose.set_psi( i, 180.0 );
		}
	}

private:
	// Private datum ///////////////////////////////////////////////////////////
	core::pose::PoseOP pose_;

};  // class RingPlaneFlipMoverTests
