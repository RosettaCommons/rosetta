// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/izstream.cxxtest.hh
/// @brief  zipstream unit test suite
/// @author Ian Davis

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

class IZStreamTests : public CxxTest::TestSuite {

	public:

	/// @brief Make sure that getline() returns the last line, even with no newline
	void test_getline() {

		std::string str;

		std::ifstream ifs( "utility/io/no_final_newline.txt" );
		int cnt_if = 0;
		//std::cout << "Reading with ifstream..." << std::endl;
		while( true ) {
			std::getline(ifs, str);
			if( !ifs ) {
				//std::cout << "good: " << ifs.good() << std::endl;
				//std::cout << "bad:  " << ifs.bad() << std::endl;
				//std::cout << "fail: " << ifs.fail() << std::endl;
				//std::cout << "eof:  " << ifs.eof() << std::endl;
				break;
			}
			//std::cout << "  " << str << std::endl;
			cnt_if++;
		}

		utility::io::izstream izs( "utility/io/no_final_newline.txt" );
		int cnt_iz = 0;
		//std::cout << "Reading with izstream..." << std::endl;
		while( true ) {
			std::getline(izs.stream(), str);
			if( !izs ) {
				//std::cout << "good: " << izs.good() << std::endl;
				//std::cout << "bad:  " << izs.bad() << std::endl;
				//std::cout << "fail: " << izs.fail() << std::endl;
				//std::cout << "eof:  " << izs.eof() << std::endl;
				break;
			}
			//std::cout << "  " << str << std::endl;
			cnt_iz++;
		}

		//std::cout << cnt_if << " lines from ifstream, " << cnt_iz << " lines from izstream" << std::endl;
		TS_ASSERT( cnt_if == cnt_iz );
	}


	void test_alternate_paths() {
		using namespace utility::io;

		izstream cannot_find( "no_final_newline.txt" );
		TS_ASSERT(!cannot_find);

		utility::vector1< std::string > current_alternative_paths(
			izstream::get_alternative_search_paths());

		utility::vector1< std::string > alternative_paths;
		alternative_paths.push_back("bla/bla/bla");
		alternative_paths.push_back("utility/io");
		alternative_paths.push_back("alb/alb/alb");

		izstream::set_alternative_search_paths(alternative_paths);

		izstream can_find( "no_final_newline.txt" );
		TS_ASSERT(can_find);

		izstream::set_alternative_search_paths(current_alternative_paths);
	}

	void test_izstream_read_from_zipped_file() {
		using namespace utility::io;
		izstream is( "utility/io/zipped_file.txt.gz" );
		std::istream & stdistream( is() );

		std::string line;
		std::getline( stdistream, line );
		//std::cout << "line: '" << line << "'" << std::endl;
		TS_ASSERT( line == "this is a file" );
		std::getline( stdistream, line );
		//std::cout << "line2: '" << line << "'" << std::endl;
		TS_ASSERT( line == "created with utmost care" );
		std::getline( stdistream, line );
		//std::cout << "line3: '" << line << "'" << std::endl;
		TS_ASSERT( line == "to be compacted" );
	}

	/// The default istream behavior on some platforms when given a directory name
	/// (rather than a filename) is to open it as a valid but empty stream.
	/// In contrast, izstream state for opening a directory should be invalid, just
	/// like it is for opening a non-existant file.
	void test_izstream_directory_read() {
		using namespace utility::io;
		// Note that for C++ streams, good() is not exactly opposite to fail() or even to bad()
		// Ain't C++ fun?

		// regular file
		izstream reg("utility/io/no_final_newline.txt");
		TS_ASSERT( reg.good() );
		TS_ASSERT( ! reg.fail() );
		// zipped file.
		izstream zip("utility/io/zipped_file.txt.gz");
		TS_ASSERT( zip.good() );
		TS_ASSERT( ! zip.fail() );
		//non-existant file.
		izstream missing("utility/io/This_file_should_not_exist._Delete_it_if_it_does.txt");
		TS_ASSERT( ! missing.good() );
		TS_ASSERT( missing.fail() );
		//directory.
		izstream dir("utility/io");
		TS_ASSERT( ! dir.good() );
		TS_ASSERT( dir.fail() );
	}

	/// izstream does not support seekg.  This is good to know.
	void dont_test_izstream_seekg() {
		using namespace utility::io;
		izstream is( "utility/io/zipped_file.txt.gz" );
		std::istream & stdistream( is() );

		std::streampos beginning( stdistream.tellg() );
		std::string line;
		std::getline( stdistream, line );
		//std::cout << "line: '" << line << "' is good? " << is.good() << std::endl;
		TS_ASSERT( line == "this is a file" );
		stdistream.seekg( beginning );
		std::getline( stdistream, line );
		//std::cout << "line, rewound: '" << line << "' is good? " << is.good() << std::endl;
		TS_ASSERT( line == "this is a file" );
		std::getline( stdistream, line );
		//std::cout << "line2: '" << line << "' is good? " << is.good()  << std::endl;
		TS_ASSERT( line == "created with utmost care" );
		std::getline( stdistream, line );
		//std::cout << "line3: '" << line << "' is good? " << is.good() << std::endl;
		TS_ASSERT( line == "to be compacted" );
	}
};

