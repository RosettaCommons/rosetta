// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/hybridization/HybridizeProtocolTests.cxxtest.hh
/// @brief  Tests for HybridizeProtocol (a Mover).  At this time there is only a trivial test on a getter/setter.
/// @author Steven Lewis (smlewi@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/hybridization/HybridizeProtocol.hh>

// Core Headers

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("HybridizeProtocolTests");


class HybridizeProtocolTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {

	}

	void tearDown() {

	}

	void test_disulf_file_getter_setter() {

		std::string const dummy_file("dummy.txt");

		core_init_with_additional_options("-fix_disulf " + dummy_file);

		using namespace protocols::hybridization;

		//The default ctor, used with XML, is supposed to read from the flag settings
		HybridizeProtocol const default_ctor;

		TS_ASSERT_EQUALS(default_ctor.get_disulf_file(), dummy_file);

		//the argumented ctor does NOT read from the flags
		//in this author's opinion that is a poor design choice (especially as this
		//setter did not exist until written for this unit test)

		//long ctor format:
		/*HybridizeProtocol(
		utility::vector1 <core::pose::PoseOP> templates_in,
		utility::vector1 <core::Real> template_weights_in,
		core::scoring::ScoreFunctionOP stage1_scorefxn_in,
		core::scoring::ScoreFunctionOP stage2_scorefxn_in,
		core::scoring::ScoreFunctionOP fa_scorefxn_in,
		std::string frag3_fn,
		std::string frag9_fn,
		std::string & cen_cst_in,
		std::string & fa_cst_in*/

		std::string const dummy_string("");
		//core::scoring::ScoreFunctionOP empty_sf(nullptr);
		utility::vector1<core::Real> const template_weights_in;
		utility::vector1<core::pose::PoseOP> const templates_in;

		std::string const frag3("protocols/fold_from_loops/movers/frags.200.3mers");
		std::string const frag9("protocols/fold_from_loops/movers/frags.200.9mers");

		//best shot at mocking in C++
		HybridizeProtocol long_ctor(templates_in, template_weights_in, nullptr, nullptr, nullptr, frag3, frag9, dummy_string, dummy_string);

		//value is unset: it should be empty instead of the option system value from above.
		TS_ASSERT_EQUALS(long_ctor.get_disulf_file(), "");

		std::string const dummy_file2("crash_test_dummy.txt");
		long_ctor.set_disulf_file(dummy_file2);
		TS_ASSERT_EQUALS(long_ctor.get_disulf_file(), dummy_file2);

	}


};
