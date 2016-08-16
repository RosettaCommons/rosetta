// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/fldsgn/filters/SheetTopologyFilter.cxxtest.hh
/// @brief  test suite for protocols::fldsgn::filters::SheetTopologyFilter
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/fldsgn/filters/SheetTopologyFilter.hh>

// Project headers

// Core headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Boost headers

// C++ headers
static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.filters.SheetTopologyFilter.cxxtest" );

void set_secstruct( core::pose::Pose & pose, std::string const & ss )
{
	debug_assert( ss.size() == pose.total_residue() );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		pose.set_secstruct( i, ss[i-1] );
	}
}
// --------------- Test Class --------------- //
class SheetTopologyFilterTest : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_residue_count_no_rgstr_shift() {
		using namespace protocols::fldsgn::filters;
		using namespace protocols::fldsgn::topology;
		SheetTopologyFilter filt;
		filt.set_secstruct( "LEEEEELHHHHHHLLHHHHHHHHHHHHHLLLEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		filt.filtered_sheet_topology( "1-2.P.99;2-3.P.99" );
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/fldsgn/filters/test_sheet.pdb" );
		core::Real score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 1.0, 1e-6 );
	}

	void test_bad_topology() {
		using namespace protocols::fldsgn::filters;
		using namespace protocols::fldsgn::topology;
		SheetTopologyFilter filt;
		filt.set_secstruct( "LEEEEELHHHHHHLLHHHHHHHHHHHHHLLLEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		filt.filtered_sheet_topology( "1-2.A.99;2-3.A.99" );
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/fldsgn/filters/test_sheet.pdb" );
		core::Real score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 0.0, 1e-4 );

		// half should be right
		filt.filtered_sheet_topology( "1-2.P.99;2-3.A.99" );
		score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 0.5, 1e-4 );

		// all should be wrong
		filt.filtered_sheet_topology( "1-3.P.99;2-3.A.99" );
		score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 0.0, 1e-4 );

		// now use actual numbers, count pairings
		filt.set_secstruct( "LEEEEEEHHHHHHLLHHHHHHHHHHHHHLLEEEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		filt.filtered_sheet_topology( "1-2.P.-1;2-3.P.0" );
		score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 0.0, 1e-4 );

		filt.set_secstruct( "LEEEEELHHHHHHLLHHHHHHHHHHHHHLLLEEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		filt.filtered_sheet_topology( "1-2.P.-1;2-3.P.0" );
		score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 0.417, 1e-2 );

		filt.set_secstruct( "LEEEEELHHHHHHLLHHHHHHHHHHHHHLLLEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		score = filt.compute( input_pose );
		TR << "Score = " << score << std::endl;
		TS_ASSERT_DELTA( score, 1.0, 1e-2 );
	}


};

