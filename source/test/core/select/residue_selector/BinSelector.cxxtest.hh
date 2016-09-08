// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/BinSelector.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::BinSelector
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/BinSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::select::residue_selector;

static THREAD_LOCAL basic::Tracer TR("BinSelectorTests");

class BinSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Confirm that the binselector properly selects residues in the alpha helix bin.
	///
	void test_binselector_abin() {
		BinSelectorOP selector( new BinSelector );
		selector->set_bin_params_file_name("ABEGO");
		selector->set_bin_name("A");
		selector->initialize_and_check();
		core::pose::Pose trpcage( create_trpcage_ideal_pose() );

		ResidueSubset subset = selector->apply( trpcage );

		for ( core::Size i=1, imax=trpcage.size(); i<=imax; ++i ) {
			TR << trpcage.residue(i).name3() << i << "\texpect:";
			if ( ( i < 2 || i > 9 ) && (i!=12 && i!=13 && i!=14) ) {
				TR << "false\t";
				TS_ASSERT( !subset[i] );
			} else {
				TR << "true\t";
				TS_ASSERT( subset[i] );
			}
			TR << "is:" << (subset[i] ? "true" : "false") << std::endl;
		}
		TR.flush();
	}

};
