// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/AngleConstraint.cxxtest.hh
/// @brief  test suite for angle constraints
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/constraints/AngleConstraint.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

// Package headers
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/conversions.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class AngleConstraintTests : public CxxTest::TestSuite
{

public:
	AngleConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_angle_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::scoring::func::FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( numeric::conversions::radians( 109 ), 10 ) );

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );

		AngleConstraint ang_cst( at1, at2, at3, func );
		EnergyMap weights, emap;
		weights[ angle_constraint ] = 1.0;
		ang_cst.score( fourpts, weights, emap );
		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ angle_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ angle_constraint ],   0.0138748361634874, 1e-14 );

		std::cout.precision( before_precision );
	}

	void test_angle_derivatives() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->size() == 2 );
		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( numeric::conversions::radians( 109 ), 10 ) );
		AngleConstraint ang_cst( at1, at2, at3, func );
		EnergyMap weights, emap;
		weights[ angle_constraint ] = 1.0;
		core::scoring::func::ConformationXYZ cfunc( ubqstump->conformation() );
		ang_cst.score( cfunc, weights, emap );
		std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ angle_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ angle_constraint ],   0.03488216167816781, 1e-14 );

		ScoreFunction sfxn;
		sfxn.set_weight( angle_constraint, 1.0 );
		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		for ( Size ii = 1; ii <= ubqstump->size(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
				at1.rsd() = at2.rsd() = at3.rsd() = ii;
				at1.atomno() = ( rsd_type.bondangle( bondang ) ).key1();
				at2.atomno() = ( rsd_type.bondangle( bondang ) ).key2();
				at3.atomno() = ( rsd_type.bondangle( bondang ) ).key3();

				AngleConstraintOP ang_cst2( new AngleConstraint( at1, at2, at3, func ) );
				ubqstump->remove_constraints();
				ubqstump->add_constraint( ang_cst2 );
				AtomDerivValidator adv( *ubqstump, sfxn, movemap );
				// This call runs a numeric deriv check on all the free dofs in the system and makes sure
				// that the analytic norm matches the numeric norm to 1e-3.
				adv.simple_deriv_check( true, 1e-6 );
			}
		}
	}

	void test_angle_constraint_clone() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( numeric::conversions::radians( 109 ), 10 ) );
		AngleConstraintOP ang_cst(  new AngleConstraint( at1, at2, at3, func ));

		ConstraintOP cloned_cst = ang_cst->clone();
		AngleConstraintOP cloned_angcst = utility::pointer::dynamic_pointer_cast< AngleConstraint > ( cloned_cst );

		// ensure the dynamic cast succeeds
		TS_ASSERT( cloned_angcst );

		// Make sure that the clone isn't the same as the original -- of course, right?
		TS_ASSERT_DIFFERS( ang_cst, cloned_cst    );
		TS_ASSERT_DIFFERS( ang_cst, cloned_angcst );

		// check mutual equality; a == b and b == a
		TS_ASSERT( *ang_cst == *cloned_cst );
		TS_ASSERT( *cloned_cst == *ang_cst );

		// clone() should perform a deep copy of the internal func object, verifiable by looking
		// at the func pointers and making sure they point at different objects.
		TS_ASSERT_DIFFERS( & ang_cst->get_func(), & cloned_cst->get_func() );

	}

	void test_angle_constraint_equality_operator() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.5, 0.5 ) );

		AngleConstraintOP angle_cst1( new AngleConstraint( at1, at2, at3, func ));

		core::scoring::func::HarmonicFuncOP func2( new core::scoring::func::HarmonicFunc( 1.5, 0.75 ) );
		AngleConstraintOP angle_cst2( new AngleConstraint( at1, at2, at3, func2 ));

		AtomID at4( 4, 1), at5( 5, 1 ), at6( 6, 1 );
		AngleConstraintOP angle_cst3( new AngleConstraint( at4, at2, at3, func ));
		AngleConstraintOP angle_cst4( new AngleConstraint( at1, at5, at3, func ));
		AngleConstraintOP angle_cst5( new AngleConstraint( at1, at2, at5, func ));

		// func objects differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *angle_cst1 != *angle_cst2 );
		TS_ASSERT( *angle_cst2 != *angle_cst1 );

		// atom1s differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *angle_cst1 != *angle_cst3 );
		TS_ASSERT( *angle_cst3 != *angle_cst1 );

		// atom2s differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *angle_cst1 != *angle_cst4 );
		TS_ASSERT( *angle_cst4 != *angle_cst1 );

		// atom3s differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *angle_cst1 != *angle_cst5 );
		TS_ASSERT( *angle_cst5 != *angle_cst1 );
	}



	void test_serialize_AngleConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AngleConstraintOP instance( new AngleConstraint( AtomID( 2,3 ), AtomID( 3,4 ), AtomID( 4,5 ), some_func ) );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		AngleConstraintOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
