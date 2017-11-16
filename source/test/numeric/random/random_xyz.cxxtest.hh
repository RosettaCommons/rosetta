// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/random_xyz.hh
/// @brief  test suite for the functions that live in numeric/random/random_xyz.hh
/// @author Darwin Fu

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

//#include <core/init_util.hh>

// Package Headers
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random_xyz.hh>

#include <numeric/constants.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("numeric.random.random_xyz_cxxtest");

class random_xyzTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init(); // Needed to setup the tracer and the RNG.
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_uniform_vector_sphere(){

		numeric::Real distance = 5;

		for ( int i=0; i <= 100; i++ ) {
			numeric::xyzVector<numeric::Real> random_point = numeric::random::uniform_vector_sphere(distance);
			TR << random_point.x() << "," << random_point.y() << "," << random_point.z() << std::endl;
			TS_ASSERT_LESS_THAN_EQUALS(random_point.length(),5);
		}
	}

};

