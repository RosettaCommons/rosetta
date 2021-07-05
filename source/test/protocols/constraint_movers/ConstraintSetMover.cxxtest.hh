// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_movers/ConstraintSetMover.cxxtest.hh
/// @brief  test for ConstraintSetMover mover
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers



// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.constraint_movers.ConstraintSetMover.cxxtest.hh");

// --------------- Test Class --------------- //

class ConstraintSetMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// This file is not built. The mover does file i/o with no option to get around it,
	// and its members can't be set by setters, so this is more trouble than it's worth
	// for me!
	void test_atompair_bounded() {
	}

};
