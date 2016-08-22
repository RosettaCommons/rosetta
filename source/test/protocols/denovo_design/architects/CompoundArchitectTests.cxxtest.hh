// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/CompoundArchitectTests.cxxtest.hh
/// @brief  Test suite for CompoundArchitect, for building de novo backbones
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/HelixArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("CompoundArchitectTests");

class CompoundArchitectTests : public CxxTest::TestSuite {
	//Define Variables
	typedef protocols::denovo_design::architects::CompoundArchitect CompoundArchitect;
	typedef protocols::denovo_design::architects::HelixArchitect HelixArchitect;
	typedef protocols::denovo_design::architects::PoseArchitect PoseArchitect;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;
	typedef protocols::denovo_design::components::StructureDataCOP StructureDataCOP;
	typedef protocols::denovo_design::connection::ConnectionArchitect ConnectionArchitect;

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}


	void test_simple_compound_architect()
	{
		std::string const arch_id = "arch";
		std::string const helix_id = "helix";
		HelixArchitect helix_arch( helix_id );
		helix_arch.set_lengths( "10" );

		CompoundArchitect arch( arch_id );
		arch.add_architect( helix_arch );

		core::Real random = 0.0250;
		StructureDataCOP sd_ptr = arch.design( core::pose::Pose(), random );
		TS_ASSERT( sd_ptr );
		StructureData const & sd = *sd_ptr;
		TS_ASSERT( sd.has_segment( helix_id ) );
		TS_ASSERT_EQUALS( sd.segment( helix_id ).elem_length(), 10 );
		TS_ASSERT_EQUALS( sd.segment( helix_id ).length(), 12 );
		TS_ASSERT_EQUALS( sd.pose_length(), 12 );
		TS_ASSERT_EQUALS( sd.length(), 10 );
		TS_ASSERT_EQUALS( sd.ss(), "LHHHHHHHHHHL" );
		TS_ASSERT_EQUALS( sd.abego(), "XAAAAAAAAAAX" );
	}

	void test_compound_architect_with_connection()
	{
		std::string const h1_id = "arch.h1";
		HelixArchitect h1( h1_id );
		h1.set_lengths( "10" );

		std::string const h2_id = "arch.h2";
		HelixArchitect h2( h2_id );
		h2.set_lengths( "12" );

		std::string const staple_id = "arch.staple";
		ConnectionArchitect conn( staple_id );
		conn.set_segment1_ids( h2_id );
		conn.set_segment2_ids( h1_id );
		conn.set_motifs( "4LX", "" );

		CompoundArchitect arch( "arch" );
		arch.add_architect( h1 );
		arch.add_architect( h2 );
		arch.add_connection( conn );

		core::Real random = 0.0250;
		StructureDataCOP sd_ptr = arch.design( core::pose::Pose(), random );
		TS_ASSERT( sd_ptr );
		StructureData const & sd = *sd_ptr;
		TR << sd << std::endl;
		TS_ASSERT( sd.has_segment( staple_id ) );
		TS_ASSERT( sd.has_segment( h1_id ) );
		TS_ASSERT( sd.has_segment( h2_id ) );
		TS_ASSERT_EQUALS( sd.segment( h1_id ).lower_segment(), staple_id );
		TS_ASSERT_EQUALS( sd.segment( h1_id ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( h2_id ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( h2_id ).upper_segment(), staple_id );
		TS_ASSERT_EQUALS( sd.segment( staple_id ).lower_segment(), h2_id );
		TS_ASSERT_EQUALS( sd.segment( staple_id ).upper_segment(), h1_id );
		TS_ASSERT_EQUALS( sd.ss(), "LHHHHHHHHHHHHLLLLHHHHHHHHHHL" );
		TS_ASSERT_EQUALS( sd.abego(), "XAAAAAAAAAAAAXXXXAAAAAAAAAAX" );
	}

	void test_gerard_case() {
		// three sidechains need to be fixed relative to one another
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/components/test_fixed_trimer.pdb" );

		CompoundArchitect comp_arch( "compound" );
		PoseArchitect pose_arch( "pose" );

		std::string const h1_id = "h1";
		HelixArchitect h1( h1_id );
		h1.set_lengths( "10" );

		std::string const h2_id = "h2";
		HelixArchitect h2( h2_id );
		h2.set_lengths( "10" );

		std::string const lig1_id = "pose.L01";

		std::string const conn1_id = "staple1";
		ConnectionArchitect conn( conn1_id );
		conn.set_motifs( "1LX", "" );
		conn.set_segment1_ids( h1_id );
		conn.set_segment2_ids( lig1_id );

		std::string const conn2_id = "staple2";
		ConnectionArchitect conn2( conn2_id );
		conn2.set_motifs( "1LX", "" );
		conn2.set_segment1_ids( lig1_id );
		conn2.set_segment2_ids( h2_id );

		comp_arch.add_architect( pose_arch );
		comp_arch.add_architect( h1 );
		comp_arch.add_architect( h2 );
		comp_arch.add_connection( conn );
		comp_arch.add_connection( conn2 );

		core::Real random = 0.2501;
		StructureDataOP sd_ptr = comp_arch.design( input_pose, random );
		TS_ASSERT( sd_ptr );
		StructureData const & sd = *sd_ptr;
		TR << "Created " << sd << std::endl;

		TS_ASSERT_EQUALS( sd.pose_length(), 31 ); // 31 = 3 (L2) + 3 (L3) + 1 (L1) + 2 (conns) + 11 (H1) + 11 (H2)
		TS_ASSERT_EQUALS( sd.segment( h1_id ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( h1_id ).upper_segment(), conn1_id );
		TS_ASSERT_EQUALS( sd.segment( conn1_id ).lower_segment(), h1_id );
		TS_ASSERT_EQUALS( sd.segment( conn1_id ).upper_segment(), lig1_id );
		TS_ASSERT_EQUALS( sd.segment( lig1_id ).lower_segment(), conn1_id );
		TS_ASSERT_EQUALS( sd.segment( lig1_id ).upper_segment(), conn2_id );
		TS_ASSERT_EQUALS( sd.segment( conn2_id ).lower_segment(), lig1_id );
		TS_ASSERT_EQUALS( sd.segment( conn2_id ).upper_segment(), h2_id );
		TS_ASSERT_EQUALS( sd.segment( h2_id ).lower_segment(), conn2_id );
		TS_ASSERT_EQUALS( sd.segment( h2_id ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( lig1_id ).movable_group(), sd.segment( "pose.L02" ).movable_group() );
		TS_ASSERT_EQUALS( sd.segment( lig1_id ).movable_group(), sd.segment( "pose.L03" ).movable_group() );
	}
};



