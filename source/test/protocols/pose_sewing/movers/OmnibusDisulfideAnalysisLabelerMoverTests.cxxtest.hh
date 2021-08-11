// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMoverTests.cxxtest.hh
/// @brief  tests the OmnibusDisulfideAnalysisLabelerMover
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/pose_sewing/movers/OmnibusDisulfideAnalysisLabelerMover.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/Remarks.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("OmnibusDisulfideAnalysisLabelerMoverTests");


class OmnibusDisulfideAnalysisLabelerMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void testOmnibusDisulfideAnalysisLabelerMover() {
		core::pose::Pose test_3_bundle = create_3_bundle_pose();
		protocols::pose_sewing::movers::OmnibusDisulfideAnalysisLabelerMoverOP mover( new protocols::pose_sewing::movers::OmnibusDisulfideAnalysisLabelerMover );
		mover->apply(test_3_bundle);
		core::Size remark_count = 0;
		utility::vector1<std::string> out_strings;
		out_strings.push_back("DSSP_DISULFIDE_1_43_L_L");
		out_strings.push_back("DISULFIDIZABLE_PAIRS: 2");
		out_strings.push_back("ALL_HELIX-HELIX_ANGLE: 126.191");
		out_strings.push_back("DISULFIDIZABLE_HELIX-HELIX_ANGLE: 126.191");
		out_strings.push_back("ALL_HELIX-HELIX_ANGLE: -34.1131");
		out_strings.push_back("ALL_HELIX-HELIX_ANGLE: -24.8175");
		out_strings.push_back("DISULFIDIZABLE_HELIX-HELIX_ANGLE: -24.8175");
		for ( auto current_remark : test_3_bundle.pdb_info()->remarks() ) {
			++remark_count;
			TS_ASSERT(current_remark.value == out_strings[remark_count]);
		}

		TS_ASSERT(remark_count == 7 );

	}

	void tearDown() {

	}






};
