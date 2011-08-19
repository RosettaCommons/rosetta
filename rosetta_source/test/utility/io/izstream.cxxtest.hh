// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/izstream.cxxtest.hh
/// @brief  zipstream unit test suite
/// @author Ian Davis

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <fstream>

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

};

