// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/VirtualSelectorTests.cxxtest.hh
/// @brief  tests the BlockSelector
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/select/residue_selector/VirtualResidueSelector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("VirtualResidueSelectorTests");


class VirtualSelectorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_virtualselector() {
		core::select::residue_selector::VirtualResidueSelectorOP selector( new core::select::residue_selector::VirtualResidueSelector );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		core::select::residue_selector::ResidueSubset subset = selector->apply( trpcage );


		for ( core::Size i=1; i<=20; ++i ) {
			TS_ASSERT( !subset[i] );
			trpcage.real_to_virtual(i);
		}
		subset = selector->apply( trpcage );
		for ( core::Size i=1; i<=20; ++i ) {
			TS_ASSERT( subset[i] );
		}

		TR.flush();
	}

	void tearDown() {

	}






};
