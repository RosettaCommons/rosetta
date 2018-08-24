// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pose/ResidueIndexDescription.hh>

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
		using namespace core::pose;

		std::string loopfile( "LOOP 2 10\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );

		auto start_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].start_res() );
		TS_ASSERT( start_rid != nullptr );
		TS_ASSERT_EQUALS( start_rid->pose_index(), 2 );

		auto end_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].end_res() );
		TS_ASSERT( end_rid != nullptr );
		TS_ASSERT_EQUALS( end_rid->pose_index(), 10 );

		auto cutpoint_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].cutpoint_res() );
		TS_ASSERT( cutpoint_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint_rid->pose_index(), 0 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		auto start2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].start_res() );
		TS_ASSERT( start2_rid != nullptr );
		TS_ASSERT_EQUALS( start2_rid->pose_index(), 12 );

		auto end2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].end_res() );
		TS_ASSERT( end2_rid != nullptr );
		TS_ASSERT_EQUALS( end2_rid->pose_index(), 18 );

		auto cutpoint2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].cutpoint_res() );
		TS_ASSERT( cutpoint2_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint2_rid->pose_index(), 0 );

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
		using namespace core::pose;

		std::string loopfile( "LOOP 2 10 8\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );

		auto start_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].start_res() );
		TS_ASSERT( start_rid != nullptr );
		TS_ASSERT_EQUALS( start_rid->pose_index(), 2 );

		auto end_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].end_res() );
		TS_ASSERT( end_rid != nullptr );
		TS_ASSERT_EQUALS( end_rid->pose_index(), 10 );

		auto cutpoint_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].cutpoint_res() );
		TS_ASSERT( cutpoint_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint_rid->pose_index(), 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.0 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		auto start2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].start_res() );
		TS_ASSERT( start2_rid != nullptr );
		TS_ASSERT_EQUALS( start2_rid->pose_index(), 12 );

		auto end2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].end_res() );
		TS_ASSERT( end2_rid != nullptr );
		TS_ASSERT_EQUALS( end2_rid->pose_index(), 18 );

		auto cutpoint2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].cutpoint_res() );
		TS_ASSERT( cutpoint2_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint2_rid->pose_index(), 0 );

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
		using namespace core::pose;

		std::string loopfile( "LOOP 2 10 8 0.3\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );

		auto start_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].start_res() );
		TS_ASSERT( start_rid != nullptr );
		TS_ASSERT_EQUALS( start_rid->pose_index(), 2 );

		auto end_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].end_res() );
		TS_ASSERT( end_rid != nullptr );
		TS_ASSERT_EQUALS( end_rid->pose_index(), 10 );

		auto cutpoint_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].cutpoint_res() );
		TS_ASSERT( cutpoint_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint_rid->pose_index(), 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.3 );
		TS_ASSERT( lfd[ 1 ].extended() == false );

		auto start2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].start_res() );
		TS_ASSERT( start2_rid != nullptr );
		TS_ASSERT_EQUALS( start2_rid->pose_index(), 12 );

		auto end2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].end_res() );
		TS_ASSERT( end2_rid != nullptr );
		TS_ASSERT_EQUALS( end2_rid->pose_index(), 18 );

		auto cutpoint2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].cutpoint_res() );
		TS_ASSERT( cutpoint2_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint2_rid->pose_index(), 0 );

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
		using namespace core::pose;

		std::string loopfile( "LOOP 2 10 8 0.3 1\nLOOP 12 18\n" );
		std::istringstream loopfstream( loopfile );
		protocols::loops::LoopsFileIO reader;
		protocols::loops::LoopsFileData lfd = *reader.read_loop_file_stream( loopfstream, "bogus" );

		TS_ASSERT( lfd.size() == 2 );

		auto start_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].start_res() );
		TS_ASSERT( start_rid != nullptr );
		TS_ASSERT_EQUALS( start_rid->pose_index(), 2 );

		auto end_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].end_res() );
		TS_ASSERT( end_rid != nullptr );
		TS_ASSERT_EQUALS( end_rid->pose_index(), 10 );

		auto cutpoint_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 1 ].cutpoint_res() );
		TS_ASSERT( cutpoint_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint_rid->pose_index(), 8 );

		TS_ASSERT( lfd[ 1 ].skip_rate() == 0.3 );
		TS_ASSERT( lfd[ 1 ].extended() == true );

		auto start2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].start_res() );
		TS_ASSERT( start2_rid != nullptr );
		TS_ASSERT_EQUALS( start2_rid->pose_index(), 12 );

		auto end2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].end_res() );
		TS_ASSERT( end2_rid != nullptr );
		TS_ASSERT_EQUALS( end2_rid->pose_index(), 18 );

		auto cutpoint2_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPoseNum const >( lfd[ 2 ].cutpoint_res() );
		TS_ASSERT( cutpoint2_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint2_rid->pose_index(), 0 );

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
		using namespace core::pose;

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

		auto start_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPDB const >( lfd[ 1 ].start_res() );
		TS_ASSERT( start_rid != nullptr );
		TS_ASSERT_EQUALS( start_rid->chain(), 'A' );
		TS_ASSERT_EQUALS( start_rid->resindex(), 2 );
		TS_ASSERT_EQUALS( start_rid->insertion_code(), ' ' );

		auto end_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPDB const >( lfd[ 1 ].end_res() );
		TS_ASSERT( end_rid != nullptr );
		TS_ASSERT_EQUALS( end_rid->chain(), 'A' );
		TS_ASSERT_EQUALS( end_rid->resindex(), 10 );
		TS_ASSERT_EQUALS( end_rid->insertion_code(), ' ' );

		auto cutpoint_rid = utility::pointer::dynamic_pointer_cast< ResidueIndexDescriptionPDB const >( lfd[ 1 ].cutpoint_res() );
		TS_ASSERT( cutpoint_rid != nullptr );
		TS_ASSERT_EQUALS( cutpoint_rid->chain(), 'A' );
		TS_ASSERT_EQUALS( cutpoint_rid->resindex(), 8 );
		TS_ASSERT_EQUALS( cutpoint_rid->insertion_code(), ' ' );

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
