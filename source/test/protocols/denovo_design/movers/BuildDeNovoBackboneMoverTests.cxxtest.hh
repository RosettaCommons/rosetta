// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/movers/BuildDeNovoBackboneMoverTests.cxxtest.hh
/// @brief  Test suite for BuildDeNovoBackboneMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/HelixArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("BuildDeNovoBackboneMoverTests");

class BuildDeNovoBackboneMoverTests : public CxxTest::TestSuite {
	//Define Variables
	typedef protocols::denovo_design::architects::CompoundArchitect CompoundArchitect;
	typedef protocols::denovo_design::architects::HelixArchitect HelixArchitect;
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
		StructureData orig_sd = StructureDataFactory::get_instance()->get_from_pose( pose );
		orig_sd.set_id( sd_id );

		std::string const conn_id = "bridge";
		ConnectionArchitect conn( conn_id );
		std::string const seg1 = "H01";
		std::string const seg2 = "H02";
		conn.set_segment1_ids( seg1 );
		conn.set_segment2_ids( seg2 );
		conn.set_motifs( "4LX", "2" );
		conn.apply( orig_sd );
		TR << "Created " << orig_sd << std::endl;

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
		orig_sd.set_cutpoint( conn_id, 0 );
		std::stringstream orig_ss;
		orig_ss << orig_sd;

		std::stringstream new_ss;
		new_ss << StructureDataFactory::get_instance()->get_from_pose( pose );

		TS_ASSERT_EQUALS( new_ss.str(), orig_ss.str() );
	}

	void test_symmetry()
	{
		using core::pose::Pose;
		using protocols::simple_moves::symmetry::SetupForSymmetryMover;

		Pose pose = create_trpcage_ideal_pose();
		Pose const input_pose = pose;

		// make symmetric
		SetupForSymmetryMover setup_sym( "protocols/denovo_design/components/C3_Z.sym" );
		setup_sym.set_preserve_datacache( true );

		std::string const sd_id = "UnitTest";
		StructureData orig_sd = StructureDataFactory::get_instance()->get_from_pose( pose );
		orig_sd.set_id( sd_id );
		TR << "orig_sd=" << orig_sd << std::endl;

		std::string const seg1 = "L02";
		std::string const seg2 = "NewHelix";

		HelixArchitect h11( seg2 );
		h11.set_lengths( "11" );

		ConnectionArchitect conn( "bridge" );
		conn.set_segment1_ids( seg1 );
		conn.set_segment2_ids( seg2 );
		conn.set_motifs( "0LX", "0" );

		BuildDeNovoBackboneMover fass;
		fass.set_score_filter( protocols::fldsgn::filters::SecondaryStructureFilter() );
		CompoundArchitect arch( sd_id );
		arch.add_architect( PoseArchitect( "trpcage" ) );
		arch.add_architect( h11 );
		arch.add_connection( conn );
		fass.add_prefold_mover( setup_sym );
		fass.set_architect( arch );
		fass.set_folder( RandomTorsionPoseFolder() );
		fass.apply( pose );

		pose.dump_pdb( "test.pdb" );
	}
};
