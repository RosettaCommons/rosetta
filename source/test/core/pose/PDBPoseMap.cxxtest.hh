// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/PDBPoseMap.cxxtest.hh
/// @brief  test for PDBPoseMap class
/// @author Steven Lewis
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>

#include <core/pose/PDBPoseMap.hh>
#include <utility/vector1.hh>

#include <test/core/init_util.hh> //necessary if there is tracer output

//Auto Headers



// --------------- Test Class --------------- //

class PDBPoseMapTests : public CxxTest::TestSuite {

public:

	// typedefs
	typedef core::pose::PDBPoseMap PDBPoseMap;

	// shared initialization
	void setUp() {
		core_init(); //necessary if there is tracer output for a failure
	}

	// shared finalization
	void tearDown() {
	}

	// ------------- Helper Functions ------------- //

	void fill_map( PDBPoseMap & map ) {
		map.insert( ' ',  -1, 'C', 1 );
		map.insert( 'X', 999, ' ', 2 );
		map.insert( 'A',  10, 'A', 3 );
	}

	// --------------- Test Cases --------------- //

	void test_PDBPoseMap_insert() {
		PDBPoseMap map;
		fill_map( map );

		TS_ASSERT_EQUALS( map.size(), 3 );
	}

	void test_PDBPoseMap_find() {
		PDBPoseMap map;
		fill_map( map );

		TS_ASSERT_EQUALS( map.find( ' ',  -1, 'C' ), 1 );
		TS_ASSERT_EQUALS( map.find( 'X', 999, ' ' ), 2 );
		TS_ASSERT_EQUALS( map.find( 'A',  10, 'A' ), 3 );
		TS_ASSERT_EQUALS( map.find( 'B',   0, ' ' ), 0 ); // 0 == not found
	}

	void test_PDBPoseMap_erase() {
		PDBPoseMap map;
		fill_map( map );
		map.erase( ' ', -1, 'C' );

		TS_ASSERT_EQUALS( map.size(), 2 );
		TS_ASSERT_EQUALS( map.find( ' ', -1, 'C' ), 0 ); // 0 == not found
	}

};//end class
