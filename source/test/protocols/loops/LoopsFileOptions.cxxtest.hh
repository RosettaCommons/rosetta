// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsFileOptions.cxxtest.hh
/// @brief test suite for protocols/loops/LoopsFileOptions.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/loops/LoopsFileOptions.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// C++ headers
#include <string>

using namespace protocols::loops;

class LoopsFileOptionsTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	// @brief test default ctor
	void test_initialize_LoopsFileOptions() {
		LoopsFileOptions opts;
		TS_ASSERT( opts.prohibit_single_residue_loops() );
	}

	void test_initialize_LoopsFileOptions_from_tags() {
		std::string loopsoptioindef = "<LoopsFileOptions name=blah prohibit_single_residue_loops=1/>\n";
		std::istringstream inputloopstream( loopsoptioindef );
		utility::tag::TagCOP tag = utility::tag::Tag::create( inputloopstream );
		LoopsFileOptions opts;
		opts.parse_my_tag( tag );
		TS_ASSERT( opts.prohibit_single_residue_loops() );
	}

	void test_initialize_LoopsFileOptions_from_tags_w_false_psrl() {
		std::string loopsoptioindef = "<LoopsFileOptions name=blah prohibit_single_residue_loops=0/>\n";
		std::istringstream inputloopstream( loopsoptioindef );
		utility::tag::TagCOP tag = utility::tag::Tag::create( inputloopstream );
		LoopsFileOptions opts;
		opts.parse_my_tag( tag );
		TS_ASSERT( ! opts.prohibit_single_residue_loops() );
	}

	void test_initialize_LoopsFileOptions_from_tags_w_default_psrl() {
		std::string loopsoptioindef = "<LoopsFileOptions name=blah/>\n";
		std::istringstream inputloopstream( loopsoptioindef );
		utility::tag::TagCOP tag = utility::tag::Tag::create( inputloopstream );
		LoopsFileOptions opts;
		opts.parse_my_tag( tag );
		TS_ASSERT( opts.prohibit_single_residue_loops() );
	}

};
