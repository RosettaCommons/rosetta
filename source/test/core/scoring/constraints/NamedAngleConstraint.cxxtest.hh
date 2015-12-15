// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/NamedAngleConstraint.cxxtest.hh
/// @brief  test suite for named angle constraints
/// @author Tom Linsky (tlinsky at uw dot edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
//#include <test/util/deriv_funcs.hh>

//Unit headers
#include <core/scoring/constraints/NamedAngleConstraint.hh>

//Protocol headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Basic headers
#include <basic/Tracer.hh>

//Numeric headers
#include <numeric/conversions.hh>

#ifdef	SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class NamedAngleConstraintTests : public CxxTest::TestSuite
{

public:
	NamedAngleConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// test score() method since I overloaded it
	void test_angle_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::Pose pose = create_trpcage_ideal_pose();

		// pose.dump_pdb( "trpcage.pdb" );
		core::scoring::func::ConformationXYZ confxyz( pose.conformation() );

		core::scoring::func::HarmonicFuncOP func(
			new core::scoring::func::HarmonicFunc( numeric::conversions::radians( 109 ), 10 ) );

		// angle is 146.3 = 2.5534167 radians
		// 109 degrees = 1.90240888 radians
		// func(x) = (x-x0)/sd = 0.065100782
		NamedAtomID at1( "N", 15 ), at2( "H", 15 ), at3( "C", 12 );
		AtomID an1( 1, 15 ), an2( 5, 15 ), an3( 3, 12 );

		AngleConstraint num_ang_cst( an1, an2, an3, func );
		NamedAngleConstraint ang_cst( at1, at2, at3, func );
		EnergyMap weights, emap, emap_num;
		weights[ angle_constraint ] = 1.0;
		ang_cst.score( confxyz, weights, emap );
		num_ang_cst.score( confxyz, weights, emap_num );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.0651, 1e-4 );
		TS_ASSERT_DELTA( emap[ angle_constraint ], emap_num[ angle_constraint ], 1e-16 );

		// angle should also be 146.3 after conversion to centroid
		// with AngleConstraint, this test will fail, but it should pass here
		core::util::switch_to_residue_type_set( pose, "centroid" );
		core::scoring::func::ConformationXYZ cenconfxyz( pose.conformation() );
		EnergyMap emap2;
		ang_cst.score( cenconfxyz, weights, emap2 );
		TS_ASSERT_DELTA( emap2[ angle_constraint ], 0.0651, 1e-4 );
	}

	void test_named_angle_constraint_clone() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 15 ), at2( 5, 15 ), at3( 3, 12 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( numeric::conversions::radians( 109 ), 10 ) );
		AngleConstraintOP ang_cst(  new AngleConstraint( at1, at2, at3, func ));

		NamedAtomID nat1( "N", 15 ), nat2( "H", 15 ), nat3( "C", 12 );
		NamedAngleConstraintOP named_ang_cst(  new NamedAngleConstraint( nat1, nat2, nat3, func ));

		ConstraintOP cloned_cst = named_ang_cst->clone();
		NamedAngleConstraintOP cloned_angcst = utility::pointer::dynamic_pointer_cast< NamedAngleConstraint > ( cloned_cst );

		// ensure the dynamic cast succeeds
		TS_ASSERT( cloned_angcst );

		// Make sure that the clone isn't the same as the original -- of course, right?
		TS_ASSERT_DIFFERS( named_ang_cst, cloned_cst    );
		TS_ASSERT_DIFFERS( named_ang_cst, cloned_angcst );

		// check mutual equality; a == b and b == a
		TS_ASSERT( *named_ang_cst == *cloned_cst );
		TS_ASSERT( *cloned_cst == *named_ang_cst );

		// clone() should perform a deep copy of the internal func object, verifiable by looking
		// at the func pointers and making sure they point at different objects.
		TS_ASSERT_DIFFERS( & named_ang_cst->get_func(), & cloned_cst->get_func() );

		// make sure that a NamedAngleConstraint is different from an AngleConstraint
		TS_ASSERT(         ang_cst->same_type_as_me( *named_ang_cst ) ); // AngleConstraint can downcast a NamedAngleConstraint to its type
		TS_ASSERT( ! named_ang_cst->same_type_as_me( *      ang_cst ) ); // but NamedAngleConstriant cannot downcast an AngleConstraint to its type

	}


	void test_serialize_NamedAngleConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		NamedAtomID nat1( "N", 15 ), nat2( "H", 15 ), nat3( "C", 12 );
		ConstraintOP instance( new NamedAngleConstraint(nat1,nat2,nat3,some_func) ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a NamedAngleConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< NamedAngleConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
