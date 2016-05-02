// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/zipstream.cxxtest.hh
/// @brief  zipstream unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/io/ocstream.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <iostream>


// --- define helper functions in separate namespace
namespace zipstream_tests {

utility::io::orstream &
output_( utility::io::orstream & os ) {

	os << "\nTest Line 1" << std::endl;
	os << "Test Line 2" << std::endl;
	return os;
}

}  // namespace zipstream_tests

using namespace zipstream_tests;


class ZipStreamTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_zipstream_general() {

		TS_ASSERT( true );// otherwise, this would be empty -- which it is, now

		// APL Disable this test because 1) it's noisy, and 2) it doesn't test to make sure anything happens.
		// print test output to standard output stream (via our wrapper)
		// output_( utility::io::oc::cout );

		// print test output to a zipped file
		// Leave out this test since there's no way to check its correctness as all it does
		// is create a file.
		//utility::io::ozstream mystream( "zipstream_test_output.gz" );
		//output_( mystream );

	}

};

