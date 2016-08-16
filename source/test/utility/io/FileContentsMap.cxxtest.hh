// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/FileContentsMap.cxxtest.hh
/// @brief  Test suite for the FileContentsMap class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/io/FileContentsMap.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

class FileContentsMapTests : public CxxTest::TestSuite {

public:
	void test_FileContentsMap_get_file_contents() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 0 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( ! fcm.has_read_limit_for_file( fname ) );

		std::string file_contents = fcm.get_file_contents( fname );
		TS_ASSERT_EQUALS( file_contents, "Testing testing 123\n" );
		TS_ASSERT( fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 1 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( fcm.has_read_limit_for_file( fname ) );

		std::string file_contents2 = fcm.get_file_contents( fname );
		TS_ASSERT_EQUALS( file_contents, file_contents2 );
		TS_ASSERT( fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 2 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( fcm.has_read_limit_for_file( fname ) );

	}

	void test_FileContentsMap_get_file_contents_ref() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 0 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( ! fcm.has_read_limit_for_file( fname ) );

		std::string const & file_contents = fcm.get_file_contents_ref( fname );
		TS_ASSERT_EQUALS( file_contents, "Testing testing 123\n" );
		TS_ASSERT( fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 1 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( fcm.has_read_limit_for_file( fname ) );

		std::string const & file_contents2 = fcm.get_file_contents_ref( fname );
		TS_ASSERT_EQUALS( file_contents, file_contents2 );
		TS_ASSERT_EQUALS( &file_contents, &file_contents2 );
		TS_ASSERT( fcm.has_file_contents( fname ) );
		TS_ASSERT_EQUALS( fcm.nreads_for_file( fname ), 2 );
		TS_ASSERT_EQUALS( fcm.nread_limit_for_file( fname ), 0 );
		TS_ASSERT( fcm.has_read_limit_for_file( fname ) );

	}


	void test_FileContentsMap_refuse_file_failure() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.refuse_unexpected_files() );
		fcm.refuse_unexpected_files( true );
		TS_ASSERT( fcm.refuse_unexpected_files() );

		try {
			fcm.get_file_contents( fname );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Unexpected file-read requested: " + fname );
		}

	}

	void test_FileContentsMap_refuse_file_success() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.refuse_unexpected_files() );
		fcm.refuse_unexpected_files( true );
		TS_ASSERT( fcm.refuse_unexpected_files() );

		fcm.increment_nread_limit( fname );

		try {
			std::string file_contents = fcm.get_file_contents( fname );
			TS_ASSERT_EQUALS( file_contents, "Testing testing 123\n" );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "oops! we shouldn't have gotten an exception here\n" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_FileContentsMap_delete_contents_at_nread_limit_1_read() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.delete_contents_at_nread_limit() );
		fcm.delete_contents_at_nread_limit( true );
		TS_ASSERT( fcm.delete_contents_at_nread_limit() );

		fcm.increment_nread_limit( fname );
		TS_ASSERT( fcm.nread_limit_for_file( fname ) == 1 );
		TS_ASSERT( fcm.nreads_for_file( fname ) == 0 );

		try {
			std::string file_contents = fcm.get_file_contents( fname );
			TS_ASSERT_EQUALS( file_contents, "Testing testing 123\n" );
			TS_ASSERT( ! fcm.has_file_contents( fname ) );
			TS_ASSERT( fcm.nreads_for_file( fname ) == 1 );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "oops! we shouldn't have gotten an exception here\n" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_FileContentsMap_delete_contents_at_nread_limit_2_read() {
		std::string fname = "utility/io/simple_input_file.txt";
		utility::io::FileContentsMap fcm;
		TS_ASSERT( ! fcm.delete_contents_at_nread_limit() );
		fcm.delete_contents_at_nread_limit( true );
		TS_ASSERT( fcm.delete_contents_at_nread_limit() );

		fcm.increment_nread_limit( fname );
		fcm.increment_nread_limit( fname );

		TS_ASSERT( fcm.nread_limit_for_file( fname ) == 2 );
		TS_ASSERT( fcm.nreads_for_file( fname ) == 0 );

		try {
			std::string file_contents1 = fcm.get_file_contents( fname );
			TS_ASSERT_EQUALS( file_contents1, "Testing testing 123\n" );
			TS_ASSERT( fcm.nreads_for_file( fname ) == 1 );
			TS_ASSERT( fcm.has_file_contents( fname ) );

			std::string file_contents2 = fcm.get_file_contents( fname );
			TS_ASSERT_EQUALS( file_contents2, "Testing testing 123\n" );
			TS_ASSERT( fcm.nreads_for_file( fname ) == 2 );
			TS_ASSERT( ! fcm.has_file_contents( fname ) );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "oops! we shouldn't have gotten an exception here\n" << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}


};

