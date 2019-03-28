// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/protocols/carbohydrates/TautomerizeAnomerMover.cxxtest.hh
/// @brief   Test suite for the TautomerizeAnomerMover
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <protocols/carbohydrates/TautomerizeAnomerMover.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic header
#include <basic/Tracer.hh>


static basic::Tracer TR( "protocols.carbohydrates.TautomerizeAnomerMover.cxxtest" );


class TautomerizeAnomerMoverTests : public CxxTest::TestSuite {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core;
		using namespace pose;

		core_init_with_additional_options( "-include_sugars" );

		pose_ = PoseOP( new Pose );
		make_pose_from_saccharide_sequence( *pose_, "b-D-Glcp-(1->4)-b-D-Glcp" );
	}

	// Destruction
	void tearDown()
	{}

	// Tests ///////////////////////////////////////////////////////////////////
	// Confirm that reducing-end and only reducing-end sugars are tautomerized.
	void test_TautomerizeAnomerMover()
	{
		using namespace std;
		using namespace core::kinematics;
		using namespace core::select::residue_selector;

		protocols::carbohydrates::TautomerizeAnomerMover mover;
		ResidueIndexSelectorOP selector( new ResidueIndexSelector );

		TR << "Testing that TautomerizeAnomerMover only tautomerizes the residues that it should." << endl;

		TR << " Residue 2: do not tautomerize (upper terminus = non-reducing end)" << endl;
		selector->set_index( "2" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT( pose_->residue( 1 ).has_property( "BETA_SUGAR" ) );
		TS_ASSERT( pose_->residue( 2 ).has_property( "BETA_SUGAR" ) );

		core::PointPosition const init_O1_xyz( pose_->residue( 1 ).xyz( "O1" ) );

		TR << " Residue 1: tautomerize (lower terminus = reducing end)" << endl;
		selector->set_index( "1" );
		mover.selector( selector );
		mover.apply( *pose_ );
		TS_ASSERT( pose_->residue( 1 ).has_property( "ALPHA_SUGAR" ) );
		TS_ASSERT( pose_->residue( 2 ).has_property( "BETA_SUGAR" ) );

		// Its property changed, but did it really move any atoms?
		TS_ASSERT_DIFFERS( pose_->residue( 1 ).xyz( "O1" ), init_O1_xyz );
	}


private:
	// Private datum ///////////////////////////////////////////////////////////
	core::pose::PoseOP pose_;

};  // class TautomerizeAnomerMoverTests
