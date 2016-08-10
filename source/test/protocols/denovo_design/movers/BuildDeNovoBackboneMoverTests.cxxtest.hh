// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/movers/BuildDeNovoBackboneMoverTests.cxxtest.hh
/// @brief  Test suite for BuildDeNovoBackboneMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("BuildDeNovoBackboneMoverTests");

class BuildDeNovoBackboneMoverTests : public CxxTest::TestSuite {
	//Define Variables
	typedef protocols::denovo_design::architects::CompoundArchitect CompoundArchitect;
	typedef protocols::denovo_design::architects::PoseArchitect PoseArchitect;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataFactory StructureDataFactory;
	typedef protocols::denovo_design::components::RandomTorsionPoseFolder RandomTorsionPoseFolder;
	typedef protocols::denovo_design::connection::ConnectionArchitect ConnectionArchitect;
	typedef protocols::denovo_design::movers::BuildDeNovoBackboneMover BuildDeNovoBackboneMover;

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}


	void test_design()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );
		// attach appropriately named StructureData
		std::string const sd_id = "UnitTest";
		StructureData orig_sd = StructureDataFactory::get_instance()->get_from_pose( pose, sd_id );

		ConnectionArchitect conn( sd_id + ".bridge" );
		std::string const seg1 = "UnitTest.H01";
		std::string const seg2 = "UnitTest.H02";
		conn.set_segment1_ids( seg1 );
		conn.set_segment2_ids( seg2 );
		conn.set_motifs( "4LX", "2" );
		conn.apply( orig_sd );
		TR << "Created " << orig_sd << std::endl;

		conn.set_id( "bridge" );
		BuildDeNovoBackboneMover fass;
		protocols::fldsgn::filters::SecondaryStructureFilter ss_filt;
		ss_filt.set_use_dssp( true );
		fass.set_score_filter( ss_filt );
		CompoundArchitect arch( sd_id );
		arch.add_architect( PoseArchitect( "twohelix" ) );
		arch.add_connection( conn );
		fass.set_architect( arch );
		fass.set_folder( RandomTorsionPoseFolder() );
		fass.apply( pose );

		// here, orig_ss has a cutpoint set -- to compare to the new SD, we need to remove it
		orig_sd.set_cutpoint( "bridge", 0 );
		std::stringstream orig_ss;
		orig_ss << orig_sd;

		std::stringstream new_ss;
		new_ss << StructureDataFactory::get_instance()->get_from_pose( pose );

		TS_ASSERT_EQUALS( new_ss.str(), orig_ss.str() );
	}
};
