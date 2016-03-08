// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/ResidueRanges.cxxtest.hh
/// @brief test suite for core::select::residue_selector::ResidueRanges
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/ResidueRanges.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.ResidueRanges.cxxtest.hh" );

class ResidueRangesTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResidueRanges_parse_my_tag()
	{
		// simple test of ranges -- residues 2-10 inclusive are selected
		core::select::residue_selector::ResidueSubset subset( 20, false );
		for ( core::Size resid=2; resid<=10; ++resid ) {
				subset[ resid ] = true;
		}

		core::select::residue_selector::ResidueRanges ranges( subset );
		TS_ASSERT_EQUALS( ranges.size(), 1 );
		TS_ASSERT_EQUALS( ranges[1].start, 2 );
		TS_ASSERT_EQUALS( ranges[1].stop, 10 );

		// now, also select the last residue
		subset[ 20 ] = true;
		ranges = core::select::residue_selector::ResidueRanges( subset );
		TS_ASSERT_EQUALS( ranges.size(), 2 );
		TS_ASSERT_EQUALS( ranges[1].start, 2 );
		TS_ASSERT_EQUALS( ranges[1].stop, 10 );
		TS_ASSERT_EQUALS( ranges[2].start, 20 );
		TS_ASSERT_EQUALS( ranges[2].stop, 20 );

		// shrink the original selection
		subset[ 2 ] = false;
		ranges = core::select::residue_selector::ResidueRanges( subset );
		TS_ASSERT_EQUALS( ranges.size(), 2 );
		TS_ASSERT_EQUALS( ranges[1].start, 3 );
		TS_ASSERT_EQUALS( ranges[1].stop, 10 );
		TS_ASSERT_EQUALS( ranges[2].start, 20 );
		TS_ASSERT_EQUALS( ranges[2].stop, 20 );

		subset[ 10 ] = false;
		ranges = core::select::residue_selector::ResidueRanges( subset );
		TS_ASSERT_EQUALS( ranges.size(), 2 );
		TS_ASSERT_EQUALS( ranges[1].start, 3 );
		TS_ASSERT_EQUALS( ranges[1].stop, 9 );
		TS_ASSERT_EQUALS( ranges[2].start, 20 );
		TS_ASSERT_EQUALS( ranges[2].stop, 20 );

		subset[ 1 ] = true;
		subset[ 2 ] = true;
		subset[ 19 ] = true;
		ranges = core::select::residue_selector::ResidueRanges( subset );
		TS_ASSERT_EQUALS( ranges.size(), 2 );
		TS_ASSERT_EQUALS( ranges[1].start, 1 );
		TS_ASSERT_EQUALS( ranges[1].stop, 9 );
		TS_ASSERT_EQUALS( ranges[2].start, 19 );
		TS_ASSERT_EQUALS( ranges[2].stop, 20 );
	}

};
