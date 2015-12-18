// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/NamedAtomPairConstraint.cxxtest.hh
/// @brief  test suite for named atom pair constraints
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
//#include <test/util/deriv_funcs.hh>

//Unit headers
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>

//Protocol headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Basic headers
#include <basic/Tracer.hh>

//Numeric headers
#include <numeric/conversions.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

static basic::Tracer TR("core.scoring.constraints.NamedAtomPairConstraint.cxxtest");

using namespace core;

class NamedAtomPairConstraintTests : public CxxTest::TestSuite
{

public:
	NamedAtomPairConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// test score() method since I overloaded it
	void test_atom_pair_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::Pose pose = create_trpcage_ideal_pose();

		// pose.dump_pdb( "trpcage.pdb" );
		core::scoring::func::ConformationXYZ confxyz( pose.conformation() );

		core::scoring::func::HarmonicFuncOP func(
			new core::scoring::func::HarmonicFunc( 1.25, .5 ) );

		// angle is 146.3 = 2.5534167 radians
		// 109 degrees = 1.90240888 radians
		// func(x) = (x-x0)/sd = 0.065100782
		NamedAtomID at1( "N", 15 ), at2( "H", 15 );
		AtomID an1( 1, 15 ), an2( 5, 15 );

		AtomPairConstraint atom_pair_cst( an1, an2, func );
		NamedAtomPairConstraint named_atom_pair_cst( at1, at2, func );
		EnergyMap weights, emap, emap_num;
		weights[ atom_pair_constraint ] = 1.0;
		atom_pair_cst.score( confxyz, weights, emap );

		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 1.0 );
		named_atom_pair_cst.setup_for_scoring( confxyz, sfxn ); // resolve the named atom constraints into actual AtomIDs
		named_atom_pair_cst.score( confxyz, weights, emap_num ); // now compute the scores

		// TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 0.0651, 1e-4 );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], emap_num[ atom_pair_constraint ], 1e-16 );

		// angle should also be 146.3 after conversion to centroid
		// with AtomPairConstraint, this test will fail, but it should pass here
		core::util::switch_to_residue_type_set( pose, "centroid" );
		core::scoring::func::ConformationXYZ cenconfxyz( pose.conformation() );
		EnergyMap emap2;
		named_atom_pair_cst.setup_for_scoring( cenconfxyz, sfxn ); // resolve the named atom constraints into actual AtomIDs
		named_atom_pair_cst.score( cenconfxyz, weights, emap2 );
		// TS_ASSERT_DELTA( emap2[ atom_pair_constraint ], 0.0651, 1e-4 );
		TS_ASSERT_DELTA( emap2[ atom_pair_constraint ], emap_num[ atom_pair_constraint ], 1e-16 );
	}

	void test_named_atom_pair_constraint_clone() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 15 ), at2( 5, 15 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.25, .5 ) );
		AtomPairConstraintOP atom_pair_cst(  new AtomPairConstraint( at1, at2, func ));

		NamedAtomID nat1( "N", 15 ), nat2( "H", 15 );
		NamedAtomPairConstraintOP named_atom_pair_cst(  new NamedAtomPairConstraint( nat1, nat2, func ));

		ConstraintOP cloned_cst = named_atom_pair_cst->clone();
		NamedAtomPairConstraintOP cloned_angcst = utility::pointer::dynamic_pointer_cast< NamedAtomPairConstraint > ( cloned_cst );

		// ensure the dynamic cast succeeds
		TS_ASSERT( cloned_angcst );

		// Make sure that the clone isn't the same as the original -- of course, right?
		TS_ASSERT_DIFFERS( named_atom_pair_cst, cloned_cst    );
		TS_ASSERT_DIFFERS( named_atom_pair_cst, cloned_angcst );

		// check mutual equality; a == b and b == a
		TS_ASSERT( *named_atom_pair_cst == *cloned_cst );
		TS_ASSERT( *cloned_cst == *named_atom_pair_cst );

		// clone() should perform a deep copy of the internal func object, verifiable by looking
		// at the func pointers and making sure they point at different objects.
		TS_ASSERT_DIFFERS( & named_atom_pair_cst->get_func(), & cloned_cst->get_func() );

		// make sure that a NamedAtomPairConstraint is different from an AtomPairConstraint
		TS_ASSERT(         atom_pair_cst->same_type_as_me( *named_atom_pair_cst ) ); // AtomPairConstraint can downcast a NamedAtomPairConstraint to its type
		TS_ASSERT( ! named_atom_pair_cst->same_type_as_me( *      atom_pair_cst ) ); // but NamedAtomPairConstriant cannot downcast an AtomPairConstraint to its type

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::func::ConformationXYZ confxyz( pose.conformation() );

		EnergyMap weights, emap, emap_num;
		ScoreFunction sfxn;
		weights[ atom_pair_constraint ] = 1.0;
		sfxn.set_weight( atom_pair_constraint, 1.0 );
		atom_pair_cst->score( confxyz, weights, emap );
		named_atom_pair_cst->setup_for_scoring( confxyz, sfxn ); // resolve the named atom constraints into actual AtomIDs
		named_atom_pair_cst->score( confxyz, weights, emap_num );

		// After scoring, the two atoms between atom_pair_cst and named_atom_pair_cst should agree, yes, but ...
		for ( core::Size ii = 1; ii <= 2; ++ii ) {
			TS_ASSERT_EQUALS( atom_pair_cst->atom(ii), named_atom_pair_cst->atom(ii) );
		}
		// .. the NamedAtomPairConstraint should still not equal the AtomPairConstraint
		TS_ASSERT( *atom_pair_cst != *named_atom_pair_cst ); // <-- this comparison would fail w/o mutual calls to same_type_as_me
		TS_ASSERT( *named_atom_pair_cst != *atom_pair_cst );
	}


	void test_serialize_NamedAtomPairConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		NamedAtomID nat1( "N", 15 ), nat2( "H", 15 );
		ConstraintOP instance( new NamedAtomPairConstraint(nat1,nat2,some_func) ); // serialize this through a pointer to the base class

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

		// make sure the deserialized base class pointer points to a NamedAtomPairConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< NamedAtomPairConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
