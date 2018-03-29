// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_movers/ResidueTypeConstraintMover.cxxtest.hh
/// @brief  test for ResidueTypeConstraintMover mover
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/constraint_movers/ResidueTypeConstraintMover.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("protocols.constraint_movers.ResidueTypeConstraintMover.cxxtest.hh");

// --------------- Test Class --------------- //

class ResidueTypeConstraintMoverTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP test_dimer_pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();

		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
	}

	void tearDown() {
	}

	void test_it_adds_constraints() {
		protocols::constraint_movers::ResidueTypeConstraintMover rtcm;

		rtcm.set_AA_name3( "ALA" );
		rtcm.set_favor_bonus( 5 );

		rtcm.apply( *test_dimer_pose_ );

		TS_ASSERT( !test_dimer_pose_->constraint_set()->is_empty() );
		TS_ASSERT( test_dimer_pose_->constraint_set()->get_all_constraints().size() == test_dimer_pose_->size() );
	}

};
