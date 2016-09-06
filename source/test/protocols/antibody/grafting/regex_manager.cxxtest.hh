// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/grafting/regex_manager.cxxtest.hh
/// @brief  test for antibody regex reading code (basically, make sure I didn't do something obviously foolish)
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#include <protocols/antibody/grafting/regex_manager.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody/grafting/util.hh>

#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting.regex_manager.cxxtest");

using std::string;

class RegExManager_tests : public CxxTest::TestSuite
{
public:
	RegExManager_tests() {}

	// Shared initialization goes here.
	void setUp() { core_init(); }

	// Shared finalization goes here.
	void tearDown() {}

	// ------------------------------------------ //
	/// @brief test if CDR detection match knownw results from our DB
	void test_cdr_regex_parsing() {
		#ifdef __ANTIBODY_GRAFTING__
			protocols::antibody::grafting::RegExManager * rem( protocols::antibody::grafting::RegExManager::get_instance() );

			// Ensure we get a chunk of memory
			TS_ASSERT( rem );

			// Test the H1 pattern
			TS_ASSERT_EQUALS( rem->H1_pattern(), "C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)" );

			// Test the H3 pattern
			TS_ASSERT_EQUALS( rem->H3_pattern(), "C[A-Z]{1,33}(W)(G|A|C)[A-Z]{1,2}(Q|S|G|R)" );

			// Test the L1 pattern
			TS_ASSERT_EQUALS( rem->L1_pattern(), "C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)" );

			// Test the L3 pattern
			TS_ASSERT_EQUALS( rem->L3_pattern(), "C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)" );
		#endif // __ANTIBODY_GRAFTING__
	}
};
