// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/Constraints.cxxtest.hh
/// @brief  test suite for the Constraints object, a container for Constraint pointers
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/types.hh>
#include <core/scoring/constraints/Constraints.hh>

#include <basic/Tracer.hh>



using basic::Error;
using basic::Warning;

//static basic::Tracer TR("core.scoring.constraints.Constraints.cxxtest");

using namespace core;
using namespace core::scoring;
using namespace core::scoring::func;
using namespace core::scoring::constraints;

class ConstraintsTests : public CxxTest::TestSuite
{

public:
	ConstraintsTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_constraints_shallow_copy_in_copy_ctor() {
		// Make sure that the copy constructor performs a shallow copy of the
		// constraints that it holds
		HarmonicFuncOP hfunc( new HarmonicFunc( 5, 0.25 ));
		AtomPairConstraintOP apc( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		Constraints csts1;
		csts1.add_constraint( apc );
		TS_ASSERT_EQUALS( csts1.size(), 1 );

		Constraints csts2( csts1 ); // copy constructor

		TS_ASSERT_EQUALS( csts2.size(), 1 );
		// shallow copy means that these two must contain pointers to the same constraint
		TS_ASSERT_EQUALS( csts1.constraints()[ 1 ], csts2.constraints()[ 1 ] );
	}

	void test_constraints_shallow_copy_in_assignment_operator() {
		// Make sure that the assignment operator performs a shallow copy of the
		// constraints that it holds
		HarmonicFuncOP hfunc( new HarmonicFunc( 5, 0.25 ));
		AtomPairConstraintOP apc( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		Constraints csts1;
		csts1.add_constraint( apc );
		TS_ASSERT_EQUALS( csts1.size(), 1 );

		Constraints csts2;
		csts2 = csts1; // assignment operator

		TS_ASSERT_EQUALS( csts2.size(), 1 );
		// shallow copy means that these two must contain pointers to the same constraint
		TS_ASSERT_EQUALS( csts1.constraints()[ 1 ], csts2.constraints()[ 1 ] );
	}

	void test_constraints_shallow_copy_in_clone() {
		// Make sure that the clone function performs a shallow copy of the
		// constraints that it holds
		HarmonicFuncOP hfunc( new HarmonicFunc( 5, 0.25 ));
		AtomPairConstraintOP apc( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		Constraints csts1;
		csts1.add_constraint( apc );
		TS_ASSERT_EQUALS( csts1.size(), 1 );

		ConstraintsOP csts2( csts1.clone() ); // clone

		TS_ASSERT_EQUALS( csts2->size(), 1 );
		// shallow copy means that these two must contain pointers to the same constraint
		TS_ASSERT_EQUALS( csts1.constraints()[ 1 ], csts2->constraints()[ 1 ] );
	}

	void test_constraints_deep_copy_in_deep_clone() {
		// Make sure that the deep_clone function performs a deep copy of the
		// constraints that it holds
		HarmonicFuncOP hfunc( new HarmonicFunc( 5, 0.25 ));
		AtomPairConstraintOP apc( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		Constraints csts1;
		csts1.add_constraint( apc );
		TS_ASSERT_EQUALS( csts1.size(), 1 );

		ConstraintsOP csts2( csts1.deep_clone() ); // deep_clone

		TS_ASSERT_EQUALS( csts2->size(), 1 );
		// deep copy means that these two must contain pointers to different constraints
		TS_ASSERT( csts1.constraints()[ 1 ] != csts2->constraints()[ 1 ] );
	}

	void test_constraints_assignment_operator_removes_old_csts() {
		// make sure that when copying old constraints into new constraints that the
		// old constraints are removed.
		HarmonicFuncOP hfunc( new HarmonicFunc( 5, 0.25 ));
		AtomPairConstraintOP apc1( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		AtomPairConstraintOP apc2( new AtomPairConstraint(
			core::id::AtomID( 1, 1 ),
			core::id::AtomID( 2, 2 )
			hfunc ));
		Constraints csts1;
		csts1.add_constraint( apc1 );
		TS_ASSERT_EQUALS( csts1.size(), 1 );

		Constraints csts2;
		csts2.add_constraint( apc2 ); // the two Constraints now hold different APCs

		csts2 = csts1; // assignment operator

		TS_ASSERT_EQUALS( csts2->size(), 1 ); // the old constraint should be gone
		TS_ASSERT_EQUALS( csts1.constraints()[ 1 ], csts2->constraints()[ 1 ] );
	}

};
