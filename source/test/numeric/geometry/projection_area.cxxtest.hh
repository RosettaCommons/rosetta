// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
//
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// /// @file test/numeric/geometry/projection_area.cxxtest.hh
// /// @brief Unit tests for projection area function.
// /// @details This is a unit test for projection area to ensure it's working the way it should and that it isn't disrupted by future additions to Rosetta.
// /// @author Bargeen Turzo (turzo.1@osu.edu)

// // Test headers
#include <cxxtest/TestSuite.h>
#include <core/init_util.hh>

// Package headers
#include <numeric/geometry/projection_area.hh>

//Project headers
#include <basic/Tracer.hh>

//Utility headers
#include <utility/vector1.hh>

static basic::Tracer TR("numeric.geometry.projection_area.cxxtest");


class ProjectionAreaTests : public CxxTest::TestSuite {

public:

	ProjectionAreaTests() {};

	//Shared initialization
	void setUp() {
		core_init();
	}

	//Shared finalization
	void tearDown() {
	}

	//test score term
	void test_projection_area() {
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace utility;

		vector1< Real > xcords { 42.375 }; // Putting xcord of an atom in vector
		vector1< Real > ycords {-12.180 }; // Putting ycord of an atom in vector
		vector1< Real > eff_radii {1.91 }; // Putting effective radius of that atom in vector
		Real probe_radius = 1.0; // Probe raidus used to approximate the area of the atom
		Real ProjectionArea= projection_area(xcords,ycords,eff_radii,probe_radius);
		TR << "Beginning ProjectionAreaTests::test_projection_area() ..." << std::endl;
		TR << "Approximate 2D projection area of an atom: " << ProjectionArea << std::endl;
		TS_ASSERT_DELTA( ProjectionArea, 9, 1e-1 ); //assess difference in projected area of the atom. It should always be 9. This is fundamentally true based on our algorithm, see more on the to be published paper.
	}
};