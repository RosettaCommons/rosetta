// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_movers/CstInfoMover.cxxtest.hh
/// @brief  test for CstInfoMover mover
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/constraint_movers/CstInfoMover.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/ConstantConstraint.hh>
#include <core/scoring/func/ConstantFunc.hh>
#include <core/pose/extra_pose_info_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("protocols.constraint_movers.CstInfoMover.cxxtest.hh");

// --------------- Test Class --------------- //

class CstInfoMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// Conditions we will want to test eventually:
	// 1. If there is a cst file provided, read in THOSE constraints; if not, use pose csts
	// 2. dump_cst_file works. (Maybe alter ConstraintIO::write_constraints to dump to stream?)
	// 3. recursive_ works in the presence of multiconstraints (as desired) and ambiguousconstraints
	void test_it_adds_scores() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::pose;

		Pose pose;

		TS_ASSERT( !hasPoseExtraScore( pose, "CST_1_measure" ) );
		TS_ASSERT( !hasPoseExtraScore( pose, "CST_1_score" ) );
		protocols::constraint_movers::CstInfoMover cim;
		pose.add_constraint( ConstraintOP( new ConstantConstraint( FuncOP( new ConstantFunc( 5 ) ) ) ) );
		cim.apply( pose );
		TS_ASSERT( hasPoseExtraScore( pose, "CST_1_measure" ) );
		TS_ASSERT( hasPoseExtraScore( pose, "CST_1_score" ) );

		core::Real value = 0;
		getPoseExtraScore( pose, "CST_1_score", value );
		TS_ASSERT( value == 5 );
	}

};
