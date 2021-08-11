// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/FractionSelectorTests.cxxtest.hh
/// @brief  tests the BlockSelector
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/select/residue_selector/FractionSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/conformation/Residue.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("FractionSelectorTests");


class FractionSelectorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_fractionselector() {
		core::select::residue_selector::FractionSelectorOP selector( new core::select::residue_selector::FractionSelector );
		core::select::residue_selector::TrueResidueSelectorOP in_selector( new core::select::residue_selector::TrueResidueSelector );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		selector->set_N(1000);
		selector->set_fraction(0.5);
		selector->set_selector( in_selector);
		core::select::residue_selector::ResidueSubset subset = selector->apply( trpcage );


		for ( core::Size i=1; i<=9; ++i ) {
			TS_ASSERT( subset[i] );
		}

		selector->set_N(5);
		subset = selector->apply( trpcage );


		for ( core::Size i=1; i<=5; ++i ) {
			TS_ASSERT( subset[i] );
		}
		selector->set_first(false);
		subset = selector->apply( trpcage );

		for ( core::Size i=16; i<=20; ++i ) {
			TS_ASSERT( subset[i] );
		}

		TR.flush();
	}

	void tearDown() {

	}






};
