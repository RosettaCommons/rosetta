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

		pose.dump_pdb( "trpcage.pdb" );
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

};
