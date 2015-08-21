// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsFileIO.cxxtest.hh
/// @brief test suite for protocols/loops/LoopsFileIO.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <string>

using core::Real;
using core::Size;
using protocols::loops::Loop;
using protocols::loops::Loops;

class LoopsFileIOTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	protocols::loops::SerializedLoopList
	posenumbered_sloops_from_loopstring(
		std::string loopfile
	)
	{
		std::istringstream lstream( loopfile );
		protocols::loops::PoseNumberedLoopFileReader reader;
		protocols::loops::SerializedLoopList sloops = reader.read_pose_numbered_loops_file( lstream, "bogus" );
		return sloops;
	}

	void test_PoseNumberedLoopFileReader_one_loop() {
		std::string loopfile( "LOOP 1 4\n" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 1 );
		TS_ASSERT( sloops[ 1 ].stop  == 4 );
		TS_ASSERT( sloops[ 1 ].cut   == 0 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
	}

	void test_PoseNumberedLoopFileReader_one_loop_w_comments() {
		std::string loopfile( "#this line will be skipped\nLOOP 1 4\n" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 1 );
		TS_ASSERT( sloops[ 1 ].stop  == 4 );
		TS_ASSERT( sloops[ 1 ].cut   == 0 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
	}

	void test_PoseNumberedLoopFileReader_two_loops() {
		std::string loopfile( "LOOP 1 4\nLOOP 5 7\n" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 2 );
		TS_ASSERT( sloops[ 1 ].start == 1 );
		TS_ASSERT( sloops[ 1 ].stop  == 4 );
		TS_ASSERT( sloops[ 1 ].cut   == 0 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
		TS_ASSERT( sloops[ 2 ].start == 5 );
		TS_ASSERT( sloops[ 2 ].stop  == 7 );
		TS_ASSERT( sloops[ 2 ].cut   == 0 );
		TS_ASSERT( sloops[ 2 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 2 ].extended == false );
	}

	void test_PoseNumberedLoopFileReader_loop_w_cutpoint() {
		std::string loopfile( "LOOP 2 15 8" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 15 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
	}

	void test_PoseNumberedLoopFileReader_loop_w_skip_rate() {
		std::string loopfile( "LOOP 2 15 8 0.3" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 15 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.3 );
		TS_ASSERT( sloops[ 1 ].extended == false );
	}

	void test_PoseNumberedLoopFileReader_loop_extended() {
		std::string loopfile( "LOOP 2 15 8 0.0 1" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 15 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == true );
	}

	void test_PoseNumberedLoopFileReader_loop_extended_old_format() {
		std::string loopfile( "LOOP 2 15 8 0.0 X" );
		protocols::loops::SerializedLoopList sloops = posenumbered_sloops_from_loopstring( loopfile );
		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 15 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == true );
	}

	void test_LoopsFileIO_read_loopstream_from_pose_numbered_file() {
		std::string loopfile( "LOOP 2 10\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 2 );

		TS_ASSERT( lfd[ 1 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 10 );

		TS_ASSERT( lfd[ 1 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_index() == 0 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		TS_ASSERT( lfd[ 2 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 12 );

		TS_ASSERT( lfd[ 2 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 18 );

		TS_ASSERT( lfd[ 2 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_index() == 0 );

		TS_ASSERT( lfd[ 2 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 2 ].extended() == false );


		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		protocols::loops::SerializedLoopList sloops = lfd.resolve_as_serialized_loops( trpcage );
		TS_ASSERT( sloops.size() == 2 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 10 );
		TS_ASSERT( sloops[ 1 ].cut   == 0 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
		TS_ASSERT( sloops[ 2 ].start == 12 );
		TS_ASSERT( sloops[ 2 ].stop  == 18 );
		TS_ASSERT( sloops[ 2 ].cut   == 0 );
		TS_ASSERT( sloops[ 2 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 2 ].extended == false );
	}

	void test_LoopsFileIO_read_loopstream_from_pose_numbered_file_w_cutpoint() {
		std::string loopfile( "LOOP 2 10 8\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 2 );

		TS_ASSERT( lfd[ 1 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 10 );

		TS_ASSERT( lfd[ 1 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_index() == 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		TS_ASSERT( lfd[ 2 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 12 );

		TS_ASSERT( lfd[ 2 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 18 );

		TS_ASSERT( lfd[ 2 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_index() == 0 );

		TS_ASSERT( lfd[ 2 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 2 ].extended() == false );


		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		protocols::loops::SerializedLoopList sloops = lfd.resolve_as_serialized_loops( trpcage );
		TS_ASSERT( sloops.size() == 2 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 10 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 1 ].extended == false );
		TS_ASSERT( sloops[ 2 ].start == 12 );
		TS_ASSERT( sloops[ 2 ].stop  == 18 );
		TS_ASSERT( sloops[ 2 ].cut   == 0 );
		TS_ASSERT( sloops[ 2 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 2 ].extended == false );
	}

	void test_LoopsFileIO_read_loopstream_from_pose_numbered_file_w_skiprate() {
		std::string loopfile( "LOOP 2 10 8 0.3\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 2 );

		TS_ASSERT( lfd[ 1 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 10 );

		TS_ASSERT( lfd[ 1 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_index() == 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.3 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		TS_ASSERT( lfd[ 2 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 12 );

		TS_ASSERT( lfd[ 2 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 18 );

		TS_ASSERT( lfd[ 2 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_index() == 0 );

		TS_ASSERT( lfd[ 2 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 2 ].extended() == false );


		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		protocols::loops::SerializedLoopList sloops = lfd.resolve_as_serialized_loops( trpcage );
		TS_ASSERT( sloops.size() == 2 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 10 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.3 );
		TS_ASSERT( sloops[ 1 ].extended == false );
		TS_ASSERT( sloops[ 2 ].start == 12 );
		TS_ASSERT( sloops[ 2 ].stop  == 18 );
		TS_ASSERT( sloops[ 2 ].cut   == 0 );
		TS_ASSERT( sloops[ 2 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 2 ].extended == false );
	}
	void test_LoopsFileIO_read_loopstream_from_pose_numbered_file_w_extended() {
		std::string loopfile( "LOOP 2 10 8 0.3 1\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 2 );

		TS_ASSERT( lfd[ 1 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 10 );

		TS_ASSERT( lfd[ 1 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_index() == 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.3 );
		TS_ASSERT( lfd[ 1 ].extended() == true );

		TS_ASSERT( lfd[ 2 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].start_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 12 );

		TS_ASSERT( lfd[ 2 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].end_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 18 );

		TS_ASSERT( lfd[ 2 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_numbered() == true );
		TS_ASSERT( lfd[ 2 ].cutpoint_res().pose_index() == 0 );

		TS_ASSERT( lfd[ 2 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 2 ].extended() == false );


		/// While we're at it, lets convert these loop indices into pose indices.

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		protocols::loops::SerializedLoopList sloops = lfd.resolve_as_serialized_loops( trpcage );
		TS_ASSERT( sloops.size() == 2 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 10 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.3 );
		TS_ASSERT( sloops[ 1 ].extended == true );
		TS_ASSERT( sloops[ 2 ].start == 12 );
		TS_ASSERT( sloops[ 2 ].stop  == 18 );
		TS_ASSERT( sloops[ 2 ].cut   == 0 );
		TS_ASSERT( sloops[ 2 ].skip_rate == 0.0 );
		TS_ASSERT( sloops[ 2 ].extended == false );
	}

	void test_LoopsFileIO_read_JSON_loops_file() {
		std::string json_loopfile =
			"# FORMAT JSON\n"
			"{\"LoopSet\" : [{\n"
			"    \"start\"  : { \"resSeq\" : 2,  \"iCode\" : \" \", \"chainID\" : \"A\" },\n"
			"    \"stop\"   : { \"resSeq\" : 10, \"iCode\" : \" \", \"chainID\" : \"A\" },\n"
			"    \"cut\"    : { \"resSeq\" : 8,  \"iCode\" : \" \", \"chainID\" : \"A\" },\n"
			"    \"extras\" : { \"skip_rate\" : 0.5, \"extend\" : true },\n"
			"  }]\n"
			"}\n";
		std::istringstream loopfstream( json_loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );
		TS_ASSERT( lfd.size() == 1 );
		TS_ASSERT( lfd[ 1 ].start_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].start_res().pose_numbered() == false );
		TS_ASSERT( lfd[ 1 ].start_res().chain() == 'A' );
		TS_ASSERT( lfd[ 1 ].start_res().resindex() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().insertion_code() == ' ' );

		TS_ASSERT( lfd[ 1 ].end_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].end_res().pose_numbered() == false );
		TS_ASSERT( lfd[ 1 ].end_res().chain() == 'A' );
		TS_ASSERT( lfd[ 1 ].end_res().resindex() == 10 );
		TS_ASSERT( lfd[ 1 ].end_res().insertion_code() == ' ' );

		TS_ASSERT( lfd[ 1 ].cutpoint_res().unassigned() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().pose_numbered() == false );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().chain() == 'A' );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().resindex() == 8 );
		TS_ASSERT( lfd[ 1 ].cutpoint_res().insertion_code() == ' ' );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.5 );
		TS_ASSERT( lfd[ 1 ].extended() == true );


		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		protocols::loops::SerializedLoopList sloops = lfd.resolve_as_serialized_loops( trpcage );

		TS_ASSERT( sloops.size() == 1 );
		TS_ASSERT( sloops[ 1 ].start == 2 );
		TS_ASSERT( sloops[ 1 ].stop  == 10 );
		TS_ASSERT( sloops[ 1 ].cut   == 8 );
		TS_ASSERT( sloops[ 1 ].skip_rate == 0.5 );
		TS_ASSERT( sloops[ 1 ].extended == true );

	}


};
