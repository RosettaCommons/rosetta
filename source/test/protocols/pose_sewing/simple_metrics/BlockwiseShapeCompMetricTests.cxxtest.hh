// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/simple_metrics/BlockwiseShapeCompMetricTests.cxxtest.hh
/// @brief  tests the BlockwiseShapeCompMetric
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers

#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <protocols/pose_sewing/simple_metrics/BlockwiseShapeCompMetric.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("BlockwiseShapeCompMetricTests");


class BlockwiseShapeCompMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_BWSCM() {
		core::pose::Pose test_3_bundle = create_3_bundle_pose();
		core::select::residue_selector::SecondaryStructureSelectorOP in_selector( new core::select::residue_selector::SecondaryStructureSelector );
		in_selector->set_selected_ss("H");
		in_selector->set_use_dssp(true);
		protocols::pose_sewing::simple_metrics::BlockwiseShapeCompMetricOP metric( new protocols::pose_sewing::simple_metrics::BlockwiseShapeCompMetric );
		metric->set_selector(in_selector);
		metric-> set_filtered_sc(0.5);
		metric-> set_filtered_area(60);
		for ( auto current_pair : metric->calculate(test_3_bundle) ) {
			TS_ASSERT(current_pair.second >= 2);
		}

	}

	void tearDown() {

	}






};
