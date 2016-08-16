// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/ConstraintSet.cxxtest.hh
/// @brief  test suite for ensuring the copy semantics of the Constraints and ConstraintSet classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// #include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/types.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/EnergyMap.hh>

// basic headers
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/conversions.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.ConstraintSet.cxxtest.hh");

using namespace core;

class ConstraintSetTests : public CxxTest::TestSuite
{

public:
	ConstraintSetTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_constraints_shallow_copy_in_assignment_op() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		ConstraintOP apc( new AtomPairConstraint( AtomID( 1, 10 ), AtomID( 2, 20 ), some_func ));
		Constraints cs1;
		cs1.add_constraint( apc );
		Constraints cs2;
		cs2 = cs1;
		ConstraintCOPs const & cs2s_csts( cs2.constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 1 );
		TS_ASSERT_EQUALS( cs2s_csts[1], apc );
	}

	void test_constraints_shallow_copy_in_copy_ctor() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		ConstraintOP apc( new AtomPairConstraint( AtomID( 1, 10 ), AtomID( 2, 20 ), some_func ));
		Constraints cs1;
		cs1.add_constraint( apc );
		Constraints cs2( cs1 );
		ConstraintCOPs const & cs2s_csts( cs2.constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 1 );
		TS_ASSERT_EQUALS( cs2s_csts[1], apc );
	}

	void test_constraints_deep_clone() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		ConstraintOP apc( new AtomPairConstraint( AtomID( 1, 10 ), AtomID( 2, 20 ), some_func ));
		Constraints cs1;
		cs1.add_constraint( apc );
		ConstraintsOP cs2 = cs1.deep_clone();
		ConstraintCOPs const & cs2s_csts( cs2->constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 1 );
		TS_ASSERT( cs2s_csts[1] != apc );
		TS_ASSERT( *cs2s_csts[1] == *apc );
	}

	void test_constraintset_assignment_operator_1() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		// 1st scenario; src CstSet has an extra constraint to a different residue that appears
		// after the constraints in the dst CstSet
		ConstraintOP apc1( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintSet cs1;
		cs1.add_constraint( apc1 );
		cs1.add_constraint( apc2 );

		ConstraintSet cs2;
		cs2.add_constraint( apc1 );
		// cs2 does not get apc2

		// copy the constraints in cs1 into cs2
		cs2 = cs1;

		ConstraintCOPs const & cs2s_csts( cs2.get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 2 );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc1 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2 ), cs2s_csts.end() );
	}

	void test_constraintset_assignment_operator_2() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		// 2nd scenario; src CstSet has constraints to a residue in the middle of the set of constraints
		// in that they both have
		ConstraintOP apc1( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintOP apc3( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 4 ), some_func ));
		ConstraintSet cs1;
		cs1.add_constraint( apc1 );
		cs1.add_constraint( apc2 );
		cs1.add_constraint( apc3 );

		ConstraintSet cs2;
		cs2.add_constraint( apc1 );
		// cs2 does not get apc2
		cs2.add_constraint( apc3 );

		// copy the constraints in cs1 into cs2
		cs2 = cs1;

		ConstraintCOPs const & cs2s_csts( cs2.get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 3 );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc1 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc3 ), cs2s_csts.end() );
	}

	void test_constraintset_assignment_operator_3() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		// 3rd scenario; dst CstSet has constraints to a residues after the end of the list of
		// csts that they both have
		ConstraintOP apc1( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintOP apc3( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 4 ), some_func ));
		ConstraintSet cs1;
		cs1.add_constraint( apc1 );
		cs1.add_constraint( apc2 );

		ConstraintSet cs2;
		cs2.add_constraint( apc1 );
		cs2.add_constraint( apc2 );
		cs2.add_constraint( apc3 );

		// copy the constraints in cs1 into cs2
		cs2 = cs1;

		ConstraintCOPs const & cs2s_csts( cs2.get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 2 );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc1 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2 ), cs2s_csts.end() );
	}

	void test_constraintset_assignment_operator_4() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		// 4th scenario; dst CstSet is missing the first cst that the src CstSet has
		ConstraintOP apc1( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintOP apc3( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 4 ), some_func ));
		ConstraintSet cs1;
		cs1.add_constraint( apc1 );
		cs1.add_constraint( apc2 );
		cs1.add_constraint( apc3 );

		ConstraintSet cs2;
		cs2.add_constraint( apc2 );
		cs2.add_constraint( apc3 );

		// copy the constraints in cs1 into cs2
		cs2 = cs1;

		ConstraintCOPs const & cs2s_csts( cs2.get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 3 );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc1 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2 ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc3 ), cs2s_csts.end() );
	}

	void test_constraintset_assignment_operator_5() {
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;
		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		// 5th scenario; dst CstSet is missing the one of the constraints between a pair of residues that
		// they both have constraints between
		ConstraintOP apc1(  new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2a( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintOP apc2b( new AtomPairConstraint( AtomID( 4, 1 ), AtomID( 4, 3 ), some_func ));
		ConstraintOP apc3(  new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 4 ), some_func ));
		ConstraintSet cs1;
		cs1.add_constraint( apc1 );
		cs1.add_constraint( apc2a );
		cs1.add_constraint( apc2b );
		cs1.add_constraint( apc3 );

		ConstraintSet cs2;
		cs2.add_constraint( apc2a );
		cs2.add_constraint( apc3 );

		// copy the constraints in cs1 into cs2
		cs2 = cs1;

		ConstraintCOPs const & cs2s_csts( cs2.get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 4 );

		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc1  ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2a ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc2b ), cs2s_csts.end() );
		TS_ASSERT_DIFFERS( std::find( cs2s_csts.begin(), cs2s_csts.end(), apc3  ), cs2s_csts.end() );
	}

	// Create a second constraint set by serializing a first one and then deserializing it; the
	// constraints contained within this constraint set should be equivalent to those in the
	// original but live in different addresses in memory.
	void test_serialize_constraint_set() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));

		ConstraintOP apc1( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 2 ), some_func ));
		ConstraintOP apc2( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 3 ), some_func ));
		ConstraintOP apc3( new AtomPairConstraint( AtomID( 4, 1 ), AtomID( 4, 3 ), some_func ));
		ConstraintOP apc4( new AtomPairConstraint( AtomID( 3, 1 ), AtomID( 3, 4 ), some_func ));
		ConstraintSetOP cs1( new ConstraintSet );

		cs1->add_constraint( apc1 );
		cs1->add_constraint( apc2 );
		cs1->add_constraint( apc3 );
		cs1->add_constraint( apc4 );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( cs1 );
		}

		std::istringstream iss( oss.str() );
		ConstraintSetOP cs2;
		{
			cereal::BinaryInputArchive arc( iss );
			arc( cs2 );
		}

		ConstraintCOPs const & cs2s_csts( cs2->get_all_constraints() );
		TS_ASSERT_EQUALS( cs2s_csts.size(), 4 );
		if ( cs2s_csts.size() != 4 ) return;

		// The constraints should be functionally identical, and
		// yet not point to the same locations in memory
		TS_ASSERT(        *apc1 == *cs2s_csts[1] );
		TS_ASSERT_DIFFERS( apc1,   cs2s_csts[1] );

		TS_ASSERT(        *apc2 == *cs2s_csts[2] );
		TS_ASSERT_DIFFERS( apc2,    cs2s_csts[2] );

		TS_ASSERT(        *apc3 == *cs2s_csts[3] );
		TS_ASSERT_DIFFERS( apc3,    cs2s_csts[3] );

		TS_ASSERT(        *apc4 == *cs2s_csts[4] );
		TS_ASSERT_DIFFERS( apc4,    cs2s_csts[4] );

#endif
	}


};
