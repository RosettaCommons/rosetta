// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/io/pdb/file_data.hh>
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

		filt.set_use_dssp( false );
		set_secstruct( input_pose, "LEEEEELHHHHHHLLHHHHHHHHHHHHHLLLEEEEELLLLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL" );
		TR << "SS=" << input_pose.secstruct() << std::endl;
		score = filt.compute( input_pose );
		TR << "1bad score = " << score << std::endl;
	}

};

