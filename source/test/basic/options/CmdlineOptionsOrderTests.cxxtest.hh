// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/options/CmdlineOptionsOrderTests.cxxtest.hh
/// @brief tests to ensure that the order of options controls option priority appropriately
/// @author Steven Lewis (smlewi@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <core/init_util.hh>

// Unit Headers
#include <core/types.hh>

// Project Headers
#include <utility/vector1.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>

// C++ Headers

//Auto Headers

static basic::Tracer TR("basic.options.CmdlineOptionsOrder");

using namespace basic::options;
using namespace basic::options::OptionKeys::fold_cst;
//convenience sample options
//-fold_cst::seq_sep_stages is a RealVector
//-fold_cst::violation_skip_basis is an Integer


/// @details the purpose of this test is to ensure that Rosetta treats repeated specifications of the same command line option properly.  vector options should append.  single-valued options take the last-specified value.  When mixed command line and flag file options are present, it behaves as if the file were expanded "in place" top to bottom.
class CmdlineOptionsOrderTests : public CxxTest::TestSuite {

public:

	void setUp() {}

	void tearDown() {}

	/// @brief does the actual TS_??? comparisons.  Same across all tests; tests vary by what setup options are specified via.
	//vsb and sss are abbreviations of the option names
	void compare_options(core::Size const vsb_ref, utility::vector1<core::Real> const & sss_ref) {

		TS_ASSERT_EQUALS(option[ violation_skip_basis ].value(), vsb_ref); //8 is passed later

		utility::vector1<core::Real> const sss_opt(option[seq_sep_stages].value());

		TS_ASSERT_EQUALS(sss_ref.size(), sss_opt.size());
		for ( core::Size i(1); i<=sss_ref.size(); ++i ) {
			TS_ASSERT_DELTA(sss_ref[i], sss_opt[i], 0.001);
		}

	}

	void test_options_order_just_cmdline() {

		utility::vector1<core::Real> const stages_ref({ 98.31, 65.75 });

		core_init_with_additional_options("-fold_cst::violation_skip_basis 3 -fold_cst::seq_sep_stages 98.31 -fold_cst::violation_skip_basis 8 -fold_cst::seq_sep_stages 65.75");

		compare_options(8, stages_ref);

	}

	void test_options_order_just_one_file() {

		utility::vector1<core::Real> const stages_ref({ 98.31, 65.75 });

		core_init_with_additional_options("@ basic/options/options_1");

		compare_options(8, stages_ref);

	}

	void test_options_order_just_flagsfiles() {
		utility::vector1<core::Real> const stages_ref({ 98.31, 65.75, 65.87, 73.19 });

		core_init_with_additional_options("@ basic/options/options_1 @ basic/options/options_2");

		compare_options(28, stages_ref);

	}

	void test_options_order_mixed() {
		utility::vector1<core::Real> const stages_ref({ 8.34, 98.31, 65.75, 52.55, 65.87, 73.19, 7777.55 });

		core_init_with_additional_options("-fold_cst::violation_skip_basis 3 -fold_cst::seq_sep_stages 8.34 @ basic/options/options_1 -fold_cst::violation_skip_basis 8 -fold_cst::seq_sep_stages 52.55 @ basic/options/options_2 -fold_cst::violation_skip_basis 888 -fold_cst::seq_sep_stages 7777.55");

		compare_options(888, stages_ref);

	}



};
