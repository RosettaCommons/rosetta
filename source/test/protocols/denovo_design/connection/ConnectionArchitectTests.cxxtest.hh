// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/ConnectionArchitectTests.cxxtest.hh
/// @brief  Unit test suite for ConnectionArchitect
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("ConnectionArchitectTests");

using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::components;
using namespace protocols::denovo_design::connection;

class ConnectionArchitectTests : public CxxTest::TestSuite {
	//Define Variables
public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_bridge_connection()
	{
		// read input pose
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_bundle_disconnected.pdb" );
		// length 15 helix
		// length 16 helix

		// setup initial structure data with "UnitTest" name
		StructureData sd = StructureDataFactory::get_instance()->get_from_pose( input_pose );

		std::string const conn_id = "bridge";
		std::string const seg1 = "H01";
		std::string const seg2 = "H02";

		ConnectionArchitect conn( conn_id );
		conn.set_segment1_ids( seg1 );
		conn.set_segment2_ids( seg2 );
		conn.set_motifs( "1LX-3EB-2LG-3EB-1LX,4LX,2LX-10HA-1LB-1LA", "3" );

		core::Real random_zero = 0.000;
		TS_ASSERT( !sd.has_segment( conn_id ) );
		conn.apply( sd, random_zero );
		TS_ASSERT( sd.has_segment( conn_id ) );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).elem_length(), 10 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).length(), 10 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).lower(), 15 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).upper(), 24 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).start(), 15 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).stop(), 24 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).cutpoint(), 17 ); // 2-->1, 3-->2, 4-->3
		TS_ASSERT_EQUALS( sd.segment( seg1 ).upper_segment(), conn_id );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).lower_segment(), seg1 );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).upper_segment(), seg2 );
		TS_ASSERT_EQUALS( sd.segment( seg2 ).lower_segment(), conn_id );
		TS_ASSERT_EQUALS( sd.segment( seg1 ).cterm_included(), true );
		TS_ASSERT_EQUALS( sd.segment( seg2 ).nterm_included(), true );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).nterm_included(), true );
		TS_ASSERT_EQUALS( sd.segment( conn_id ).cterm_included(), true );

		// this should fail since the groups are the same
		conn.set_motifs( "1LX-3EB-2LG-3EB-1LX,4LX,2LX-10HA-1LB-1LA", "" );
		sd = StructureDataFactory::get_instance()->get_from_pose( input_pose );
		TS_ASSERT_THROWS( conn.apply( sd, random_zero ), EXCN_ConnectionSetupFailed );
	}

	void test_staple_connection()
	{
		// read input pose
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_bundle_disconnected.pdb" );

		// set movable groups so they are different
		StructureData orig_sd = StructureDataFactory::get_instance()->get_from_pose( input_pose );
		std::string const seg1 = "H01";
		std::string const seg2 = "H02";
		TS_ASSERT( orig_sd.has_segment( seg2 ) );
		orig_sd.set_movable_group( seg2, 2 );
		StructureDataFactory::get_instance()->save_into_pose( input_pose, orig_sd );

		StructureData pose_sd = StructureDataFactory::get_instance()->get_from_pose( input_pose );
		TS_ASSERT_EQUALS( pose_sd.segment( seg1 ).movable_group(), 1 );
		TS_ASSERT_EQUALS( pose_sd.segment( seg2 ).movable_group(), 2 );

		std::string const conn_id = "staple";
		ConnectionArchitect conn( conn_id );
		conn.set_segment1_ids( seg1 );
		conn.set_segment2_ids( seg2 );
		conn.set_motifs( "1LX-3EB-2LG-3EB-1LX,4LX,2LX-10HA-1LB-1LA", "" );

		core::Real random_zero = 0.000;
		conn.apply( pose_sd, random_zero );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).elem_length(), 10 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).length(), 10 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).lower(), 15 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).upper(), 24 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).start(), 15 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).stop(), 24 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).cutpoint(), 0 ); // 2-->1, 3-->2, 4-->3
		TS_ASSERT_EQUALS( pose_sd.segment( seg1 ).upper_segment(), conn_id );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).lower_segment(), seg1 );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).upper_segment(), seg2 );
		TS_ASSERT_EQUALS( pose_sd.segment( seg2 ).lower_segment(), conn_id );
		TS_ASSERT_EQUALS( pose_sd.segment( seg1 ).cterm_included(), true );
		TS_ASSERT_EQUALS( pose_sd.segment( seg2 ).nterm_included(), true );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).nterm_included(), true );
		TS_ASSERT_EQUALS( pose_sd.segment( conn_id ).cterm_included(), true );
	}
};

