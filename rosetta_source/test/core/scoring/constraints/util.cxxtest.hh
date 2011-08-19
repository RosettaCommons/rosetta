// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/protein_silent.cxxtest.hh
/// @brief  test suite for protein silent-file format
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/scoring/constraints/util.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <utility/stream_util.hh>


using namespace core;

class ConstraintUtilTests : public CxxTest::TestSuite {
public:
	ConstraintUtilTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_gaussian_functions() {
		using core::Size;
	 	using core::Real;
		using utility::vector1;
		using namespace core::scoring::constraints;

		Real const TOLERATED_ERROR( 1e-2 );

		Real const mean( 0 );
		Real const sdev( 1 );
		Real const weight( 1 );

		TS_ASSERT_DELTA(
			dgaussian( 0, mean, sdev, weight ), 0.3969, TOLERATED_ERROR
		);
		TS_ASSERT_DELTA(
			logdgaussian( 0, mean, sdev, weight ), -0.9239, TOLERATED_ERROR
		);

	}

	void test_exponential_functions() {
		using core::Size;
	 	using core::Real;
		using utility::vector1;
		using namespace core::scoring::constraints;

		Real const TOLERATED_ERROR( 1e-2 );

		Real const anchor( 0 );
		Real const rate  ( 1 );
		Real const weight( 1 );

		TS_ASSERT_DELTA(
			dexponential( 1, anchor, rate, weight ), 0.3678, TOLERATED_ERROR
		);
	}
}; // ConstraintUtilTestsTests
