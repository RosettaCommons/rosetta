// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/sewing/SewingEdgeFileReaderTests.cxxtest.hh
/// @brief  a unit test to ensure proper edge file importing
/// @author Minnie Langlois (minnie@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/sewing/sampling/EdgeFileReader.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>


using namespace protocols::sewing;
static basic::Tracer TR("SewingEdgeFileReaderTests");


class SewingEdgeFileReaderTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){
	}

	void test_data_members_arent_set_by_default_constructor(){
		EdgeFileReader edge_file_reader = EdgeFileReader();
		TS_ASSERT_EQUALS(edge_file_reader.is_edge_file_loaded(), false);
	}

	void test_bad_edge_file_throws_error(){
		std::string bad_edge_file = "thisfiledoesnotexist.txt";
		try{
			EdgeFileReader edge_file_reader = EdgeFileReader();
			edge_file_reader.load_edge_file( bad_edge_file );
			TS_FAIL( "Bad edge file doesn't throw error!" );
		} catch (...) {}
	}

	void test_edge_file_reader_imports_edge_map(){
		EdgeMap dummy_edge_map;

		core::Size seg_1 = 1;
		core::Size seg_3 = 3;

		std::set< std::pair< int, core::Size > > seg11_set;
		seg11_set.insert(std::make_pair( 2, seg_3 ) );
		dummy_edge_map[ std::make_pair( 1, seg_1) ] = seg11_set;

		std::set< std::pair< int, core::Size > > seg13_set;
		seg13_set.insert(std::make_pair( 3, seg_1 ) );
		dummy_edge_map[ std::make_pair( 1, seg_3) ] = seg13_set;

		std::set< std::pair< int, core::Size > > seg23_set;
		seg23_set.insert(std::make_pair( 1, seg_1 ) );
		dummy_edge_map[ std::make_pair( 2, seg_3 ) ] = seg23_set;

		std::set< std::pair< int, core::Size > > seg31_set;
		seg31_set.insert(std::make_pair( 1, seg_3 ) );
		seg31_set.insert(std::make_pair( 4, seg_3 ) );
		dummy_edge_map[ std::make_pair( 3, seg_1 ) ] = seg31_set;

		std::set< std::pair< int, core::Size > > seg33_set;
		seg33_set.insert(std::make_pair( 5, seg_1 ) );
		dummy_edge_map[ std::make_pair( 3, seg_3 ) ] = seg33_set;

		std::set< std::pair< int, core::Size > > seg43_set;
		seg43_set.insert(std::make_pair( 3, seg_1 ) );
		dummy_edge_map[ std::make_pair( 4, seg_3 ) ] = seg43_set;

		std::set< std::pair< int, core::Size > > seg51_set;
		seg51_set.insert(std::make_pair( 3, seg_3 ) );
		dummy_edge_map[ std::make_pair( 5, seg_1 ) ] = seg51_set;
		try {
			EdgeFileReader edge_file_reader = EdgeFileReader( dummy_edge_file_name );
			TS_ASSERT( edge_file_reader.is_edge_file_loaded() );

			// check edge map
			TS_ASSERT_EQUALS( edge_file_reader.edge_map(), dummy_edge_map );
			TS_ASSERT_EQUALS( edge_file_reader.get_edges( std::make_pair( 3, seg_1 ) )->second , seg31_set );
			TS_ASSERT_EQUALS( edge_file_reader.get_edges( std::make_pair( 4, seg_3 ) )->second , seg43_set );

			// check settings
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().max_clash_score, 2 );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().min_hash_score, 20 );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().boxes_per_dimension, 3 );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().hash_between_termini, 1 );
		} catch (...){
			TS_FAIL( "Could not load dummy edge file!" );
		}
	}

	void test_copied_edge_file_readers_are_equal(){
		try {
			EdgeFileReader edge_file_reader = EdgeFileReader( dummy_edge_file_name );
			EdgeFileReader edge_file_reader_copy = EdgeFileReader( edge_file_reader );

			TS_ASSERT_EQUALS( edge_file_reader.is_edge_file_loaded(), edge_file_reader_copy.is_edge_file_loaded() );
			TS_ASSERT_EQUALS( edge_file_reader.edge_map(), edge_file_reader_copy.edge_map() );
			TS_ASSERT_EQUALS( edge_file_reader.version(), edge_file_reader_copy.version() );
			TS_ASSERT_EQUALS( edge_file_reader.edge_file(), edge_file_reader_copy.edge_file() );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().max_clash_score, edge_file_reader_copy.hasher_settings().max_clash_score );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().min_hash_score, edge_file_reader_copy.hasher_settings().min_hash_score );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().boxes_per_dimension, edge_file_reader_copy.hasher_settings().boxes_per_dimension );
			TS_ASSERT_EQUALS( edge_file_reader.hasher_settings().hash_between_termini, edge_file_reader_copy.hasher_settings().hash_between_termini );

		} catch (...){
			TS_FAIL( "Could not load dummy edge file!" );
		}
	}

private:
	std::string dummy_edge_file_name = "protocols/sewing/inputs/dummy_edge_file";
};



