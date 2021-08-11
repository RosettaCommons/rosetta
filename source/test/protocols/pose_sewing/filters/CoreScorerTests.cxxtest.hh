// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/filters/CoreScorerTests.cxxtest.hh
/// @brief  tests the CoreScorer
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers

#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <protocols/pose_sewing/filters/CoreScorer.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("CoreScorerTests");


class CoreScorerTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_corescorer() {
		core::select::residue_selector::SecondaryStructureSelectorOP in_selector( new core::select::residue_selector::SecondaryStructureSelector );
		in_selector->set_selected_ss("H");
		in_selector->set_use_dssp(true);
		core::pose::Pose test_3_bundle = create_3_bundle_pose();
		protocols::pose_sewing::filters::CoreScorerOP filter(new protocols::pose_sewing::filters::CoreScorer);

		filter->set_selector(in_selector);
		filter->set_score_cutoff(-0.2);
		filter->set_window_width(6);
		filter->set_link_cutoff(3);
		filter->set_sum(true);
		filter->set_distance_mode(true);
		filter->set_distance_cutoff(10);

		//TS_ASSERT(filter->report_sm(test_3_bundle) > 0);
		TS_ASSERT( filter->apply(test_3_bundle) );
	}

	void tearDown() {

	}






};
