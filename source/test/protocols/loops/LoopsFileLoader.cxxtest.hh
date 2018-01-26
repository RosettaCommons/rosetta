// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileLoader.cxxtest.hh
/// @brief  test suite for protocols::loops::LoopsFileLoader
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/loops/LoopsFileLoader.hh>
#include <protocols/loops/LoopsFileIO.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// basic headers
#include <basic/resource_manager/ResourceManager.hh>

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
		LoopsFileLoader loader;
		std::string loopfile( "LOOP 1 4\nLOOP 5 7\n" );
		std::istringstream lstream1( loopfile );
		std::istringstream lstream2( loopfile );
		basic::resource_manager::ResourceManager rm;
		utility::tag::TagCOP tag1 = utility::tag::Tag::create( "<LoopsFile name=\"blah\" prohibit_single_residue_loops=\"true\"/>\n" );
		utility::tag::TagCOP tag2 = utility::tag::Tag::create( "<LoopsFile name=\"blah\" prohibit_single_residue_loops=\"false\"/>\n" );

		utility::pointer::ReferenceCountCOP resource1 = loader.create_resource( rm, tag1, "unit_test", lstream1 );
		utility::pointer::ReferenceCountCOP resource2 = loader.create_resource( rm, tag2, "unit_test", lstream2 );
		TS_ASSERT( resource1 ); // make sure a resource was returned
		TS_ASSERT( resource2 ); // make sure a resource was returned

		LoopsFileDataCOP lfdptr1 = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData const > ( resource1 );
		LoopsFileDataCOP lfdptr2 = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData const > ( resource2 );
		TS_ASSERT( lfdptr1 ); // make sure we're actually returned the correct type
		TS_ASSERT( lfdptr2 ); // make sure we're actually returned the correct type

		{
			LoopsFileData const & lfd( *lfdptr1 );
			TS_ASSERT( lfd.size() == 2 );
			TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 1 );
			TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 4 );
			TS_ASSERT( lfd[ 1 ].prohibit_single_residue_loops() );
			TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 5 );
			TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 7 );
			TS_ASSERT( lfd[ 2 ].prohibit_single_residue_loops() );
		}

		{
			LoopsFileData const & lfd( *lfdptr2 );
			TS_ASSERT( lfd.size() == 2 );
			TS_ASSERT( lfd[ 1 ].start_res().pose_index() == 1 );
			TS_ASSERT( lfd[ 1 ].end_res().pose_index() == 4 );
			TS_ASSERT( ! lfd[ 1 ].prohibit_single_residue_loops() );
			TS_ASSERT( lfd[ 2 ].start_res().pose_index() == 5 );
			TS_ASSERT( lfd[ 2 ].end_res().pose_index() == 7 );
			TS_ASSERT( ! lfd[ 2 ].prohibit_single_residue_loops() );
		}

	}


};
