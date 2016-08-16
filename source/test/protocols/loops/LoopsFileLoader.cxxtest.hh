// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopsFileOptions.cxxtest.hh
/// @brief test suite for protocols/loops/LoopsFileOptions.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/loops/LoopsFileOptions.hh>
#include <protocols/loops/LoopsFileLoader.hh>
#include <protocols/loops/LoopsFileIO.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Numeric headers

// C++ headers
#include <string>

using namespace protocols::loops;

class LoopsFileLoaderTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	// @brief test default ctor
	void test_create_LoopsFileData_from_LoopsFileLoder() {
		LoopsFileOptions opts;
		LoopsFileLoader loader;
		std::string loopfile( "LOOP 1 4\nLOOP 5 7\n" );
		std::istringstream lstream( loopfile );

		utility::pointer::ReferenceCountOP resource = loader.create_resource( opts, "unit_test", lstream );
		TS_ASSERT( resource ); // make sure a resource was returned

		LoopsFileDataOP lfdptr = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData > ( resource );
		TS_ASSERT( lfdptr ); // make sure we're actually returned the correct type

		LoopsFileData const & lfd( *lfdptr );
		TS_ASSERT( lfd.size() == 2 );
		TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 1 );
		TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 4 );
		TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 5 );
		TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 7 );
	}


};
