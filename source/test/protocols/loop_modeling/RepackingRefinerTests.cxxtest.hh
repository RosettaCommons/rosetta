// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/rosettascripts.hh>

// Protocol headers
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>

// C++ headers
#include <string>

// Namespaces
using namespace std;
using namespace protocols::loop_modeling;
using protocols::loop_modeling::refiners::RepackingRefiner;
using protocols::loop_modeling::refiners::RepackingRefinerOP;

class RepackingRefinerTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}

	void test_repacking_refiner_options() {
		string tag = "<RepackingRefiner once_every=50/>";
		RepackingRefinerOP refiner = parse_tag<RepackingRefiner>(tag);

		TS_ASSERT_EQUALS(refiner->get_repack_period(), 50);
	}


};

