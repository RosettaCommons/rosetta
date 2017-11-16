// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/fldsgn/filters/HelixPairingFilterTests.cxxtest.hh
/// @brief  Test suite for HelixPairingFilter
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <protocols/fldsgn/filters/HelixPairingFilter.hh>
#include <protocols/moves/DsspMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/build_pose_as_is.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("HelixPairingFilterTests");

using namespace protocols::fldsgn::filters;

class HelixPairingFilterTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_pose_without_helices()
	{
		// pose without helices should be handled gracefully
		core::pose::Pose pose = create_trpcage_ideal_pose();

		HelixPairingFilter filt( "1-2.P;2-3.A" );
		filt.secstruct( std::string( pose.size(), 'L' ) );

		// filter should fail, because helices 1, 2 and 3 are given, and the pose doesn't contain helices
		TS_ASSERT( !filt.apply( pose ) );
		TR << "Filter returns " << filt.apply( pose ) << std::endl;
		TR << "Filter metric " << filt.report_sm( pose ) << std::endl;

		// filter should pass, because no helix pairings are given, so there is nothing to check
		HelixPairingFilter filt2;
		TS_ASSERT( filt2.apply( pose ) );
	}

	void test_two_helix()
	{
		// Pose properties
		// dist = 7.9 A
		// bend_angle = 29
		// cross_angle = 8
		// align_angle = -99
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );
		protocols::moves::DsspMover().apply( pose );

		// This should pass
		HelixPairingFilter hpair( "1-2.A" );
		hpair.bend_angle( 30.0 );
		TS_ASSERT( hpair.apply( pose ) );

		// This should fail
		hpair.helix_pairings( "1-2.P" );
		TS_ASSERT( !hpair.apply( pose ) );
	}
};



