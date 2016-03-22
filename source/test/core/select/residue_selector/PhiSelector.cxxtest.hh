// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/PhiSelector.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::PhiSelector
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/PhiSelector.hh>

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

static THREAD_LOCAL basic::Tracer TR("PhiSelectorTests");

class PhiSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Confirm that the phiselector properly selects residues in the positive phi space.
	///
	void test_phiselector_positivephi() {
		PhiSelectorOP selector( new PhiSelector );
		core::pose::Pose trpcage( create_trpcage_ideal_pose() );

		ResidueSubset subset = selector->apply( trpcage );

		for ( core::Size i=1, imax=trpcage.n_residue(); i<=imax; ++i ) {
			TR << trpcage.residue(i).name3() << i << "\texpect:";
			if ( i!=10 && i!=11 && i!=15 ) {
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

	/// @brief Confirm that the phiselector properly selects residues in the negative phi space.
	///
	void test_phiselector_negativephi() {
		PhiSelectorOP selector( new PhiSelector );
		selector->set_select_positive_phi(false);
		core::pose::Pose trpcage( create_trpcage_ideal_pose() );

		ResidueSubset subset = selector->apply( trpcage );

		for ( core::Size i=1, imax=trpcage.n_residue(); i<=imax; ++i ) {
			TR << trpcage.residue(i).name3() << i << "\texpect:";
			if ( i!=1 && i!=10 && i!=11 && i!=15 && i!=imax ) {
				TR << "true\t";
				TS_ASSERT( subset[i] );
			} else {
				TR << "false\t";
				TS_ASSERT( !subset[i] );
			}
			TR << "is:" << (subset[i] ? "true" : "false") << std::endl;
		}
		TR.flush();
	}

	/// @brief Confirm that the phiselector properly selects residues in the negative phi space when termini can be selected.
	///
	void test_phiselector_negativephi_no_ignore_unconnected_upper() {
		PhiSelectorOP selector( new PhiSelector );
		selector->set_select_positive_phi(false);
		selector->set_ignore_unconnected_upper(false);
		core::pose::Pose trpcage( create_trpcage_ideal_pose() );

		ResidueSubset subset = selector->apply( trpcage );

		for ( core::Size i=1, imax=trpcage.n_residue(); i<=imax; ++i ) {
			TR << trpcage.residue(i).name3() << i << "\texpect:";
			if ( i!=1 && i!=10 && i!=11 && i!=15 ) {
				TR << "true\t";
				TS_ASSERT( subset[i] );
			} else {
				TR << "false\t";
				TS_ASSERT( !subset[i] );
			}
			TR << "is:" << (subset[i] ? "true" : "false") << std::endl;
		}
		TR.flush();
	}

};
