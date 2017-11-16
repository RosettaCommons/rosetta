// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/simple_moves/SmallMoverTests.cxxtest.hh
/// @brief  Unit tests for the SmallMover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //Convenience -- for building peptides
#include <protocols/simple_moves/BackboneMover.hh>
#include <numeric/angle.functions.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SmallMoverTests");

#define ERROR_MARGIN 0.00001

class SmallMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


	/// @brief Test that the SmallMover works properly with oligoureas.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_perturb_oligourea(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_VAL", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi( 2, -63.1 );
		pose.set_theta( 2, 75.7 );
		pose.set_psi( 2, -33.9 );
		pose.set_mu( 2, 178.3 );
		pose.set_omega( 2, -177.4 );

		TS_ASSERT_DELTA( pose.phi(2), -63.1, ERROR_MARGIN );
		TS_ASSERT_DELTA( pose.theta(2), 75.7, ERROR_MARGIN );
		TS_ASSERT_DELTA( pose.psi(2), -33.9, ERROR_MARGIN );
		TS_ASSERT_DELTA( pose.mu(2), 178.3, ERROR_MARGIN );
		TS_ASSERT_DELTA( pose.omega(2), -177.4, ERROR_MARGIN );

		protocols::simple_moves::SmallMover small;
		small.angle_max( 12.0 );
		small.nmoves(2000);
		small.temperature(100.0);
		small.apply(pose);

		TS_ASSERT_LESS_THAN( 0.05, std::abs( numeric::principal_angle_degrees( pose.phi(2)   + 63.1 ) ) );
		TS_ASSERT_LESS_THAN( 0.05, std::abs( numeric::principal_angle_degrees( pose.theta(2) - 75.7 ) ) );
		TS_ASSERT_LESS_THAN( 0.05, std::abs( numeric::principal_angle_degrees( pose.psi(2)   + 33.9 ) ) );
		TS_ASSERT_DELTA( pose.mu(2), 178.3, ERROR_MARGIN );
		TS_ASSERT_DELTA( pose.omega(2), -177.4, ERROR_MARGIN );
	}



};
