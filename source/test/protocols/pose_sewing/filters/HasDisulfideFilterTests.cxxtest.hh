// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/filters/HasDisulfideFilterTests.cxxtest.hh
/// @brief  tests the HasDisulfideFilter
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/pose_sewing/filters/HasDisulfideFilter.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/FalseResidueSelector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("HasDisulfideFilterTests");


class HasDisulfideFilterTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_hasDisulfideFilter() {
		core::select::residue_selector::TrueResidueSelectorOP in_selector( new core::select::residue_selector::TrueResidueSelector );
		core::select::residue_selector::FalseResidueSelectorOP not_in_selector( new core::select::residue_selector::FalseResidueSelector );
		core::pose::Pose test_3_bundle = create_3_bundle_pose();
		protocols::pose_sewing::filters::HasDisulfideFilterOP filter(new protocols::pose_sewing::filters::HasDisulfideFilter);

		filter->set_first_selector(in_selector);
		filter->set_second_selector(in_selector);

		TS_ASSERT(filter->apply(test_3_bundle));

		filter->set_second_selector(not_in_selector);

		TS_ASSERT(!filter->apply(test_3_bundle));
	}

	void tearDown() {

	}






};
