// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_filters/ConstraintScoreCutoffFilter.cxxtest.hh
/// @brief  test for ConstraintScoreCutoffFilter mover
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/constraint_filters/ConstraintScoreCutoffFilter.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/ConstantConstraint.hh>
#include <core/scoring/func/ConstantFunc.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("protocols.constraint_filters.ConstraintScoreCutoffFilter.cxxtest.hh");

// --------------- Test Class --------------- //

class ConstraintScoreCutoffFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// Conditions we will want to test eventually:
	// 1. We should be able to use this filter (but can't yet) on the pose's OWN CONSTRAINTS.
	void test_it_adds_scores() {
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::pose;

		Pose pose;

		ConstraintCOPs csts;
		csts.emplace_back( new ConstantConstraint( FuncOP( new ConstantFunc( 5 ) ) ) );
		protocols::constraint_filters::ConstraintScoreCutoffFilter cim( 5 );
		cim.set_constraints( csts );
		cim.set_score_type( constant_constraint );
		TS_ASSERT( cim.apply( pose ) );

		cim.set_cutoff( 4 );
		TS_ASSERT( !cim.apply( pose ) );

		cim.set_cutoff( 6 );
		TS_ASSERT( cim.apply( pose ) );
	}

};
