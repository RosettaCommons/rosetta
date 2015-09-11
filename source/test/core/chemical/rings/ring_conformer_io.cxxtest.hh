// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/chemical/rings/ring_conformer_io.cxxtest.hh
/// @brief   Test suite for ring conformer database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/rings/ring_conformer_io.hh>
#include <core/chemical/rings/RingConformerSet.hh>

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <map>


class RingConformerIOTests : public CxxTest::TestSuite {
public:
	// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


	// Confirm that ring conformers are loaded correctly from the database.
	void test_conformers_from_database_file_for_ring_size()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical::rings;

		TS_TRACE( "Testing read_conformers_from_database_file_for_ring_size() method." );

		vector1< RingConformer > const conformers(
				read_conformers_from_database_file_for_ring_size( "core/chemical/rings/dummy_conformers.data", 8 ) );

		TS_ASSERT_EQUALS( conformers.size(), 3);
		TS_ASSERT_EQUALS( conformers[ 1 ].specific_name, "1F2" );
		TS_ASSERT_EQUALS( conformers[ 2 ].general_name, "bar" );
		TS_ASSERT_EQUALS( conformers[ 3 ].degeneracy, 3 );
		TS_ASSERT_EQUALS( conformers[ 1 ].CP_parameters.size(), 5 );  // 3 fewer parameters than the ring size.
		TS_ASSERT_EQUALS( conformers[ 2 ].CP_parameters[ q ], 0.5 );
		TS_ASSERT_EQUALS( conformers[ 3 ].CP_parameters[ PHI ], 135.0 );
		TS_ASSERT_EQUALS( conformers[ 1 ].CP_parameters[ THETA ], 30.0 );
		TS_ASSERT_EQUALS( conformers[ 2 ].CP_parameters[ 5 ], 120.0 );
		TS_ASSERT_EQUALS( conformers[ 3 ].nu_angles.size(), 7 );  // 1 fewer angles than the ring size should be read.
		TS_ASSERT_EQUALS( conformers[ 1 ].nu_angles[ 1 ], -60.0 );
		TS_ASSERT_EQUALS( conformers[ 2 ].nu_angles[ 3 ], 60.0 );
		TS_ASSERT_EQUALS( conformers[ 3 ].nu_angles[ 5 ], 30.0 );
		TS_ASSERT_EQUALS( conformers[ 3 ].tau_angles.size(), 8 );  // same number of angles as ring size should be read.
		TS_ASSERT_EQUALS( conformers[ 1 ].tau_angles[ 1 ], 123.4 );
		TS_ASSERT_EQUALS( conformers[ 2 ].tau_angles[ 4 ], 345.6 );
		TS_ASSERT_EQUALS( conformers[ 3 ].tau_angles[ 6 ], 123.4 );
	}
};  // class RingConformerIOTests
