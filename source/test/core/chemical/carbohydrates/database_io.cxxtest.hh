// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 database_io.cxxtest.hh
/// @brief   Test suite for carbohydrate database loading
/// @author  Labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/carbohydrates/database_io.hh>
#include <core/chemical/carbohydrates/RingConformerSet.hh>

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <map>


class CarbohydrateDatabaseIOTests : public CxxTest::TestSuite {
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


	// Tests //////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Confirm that carbohydrate properties are loaded correctly from the database.
	void test_read_properties_from_database_file()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical::carbohydrates;

		TS_TRACE("Testing read_properties_from_database_file() method.");

		vector1<string> properties =
				read_properties_from_database_file("core/chemical/carbohydrates/sugar_properties.list");

		TS_ASSERT_EQUALS(properties.size(), 19);
	}

	// Confirm that carbohydrate 3-letter codes and roots are loaded correctly from the database.
	void test_read_codes_and_roots_from_database_file()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical::carbohydrates;

		TS_TRACE("Testing read_codes_and_roots_from_database_file() method.");

		map<string, string> map =
				read_codes_and_roots_from_database_file("core/chemical/carbohydrates/codes_to_roots.map");

		TS_ASSERT_EQUALS(map.size(), 3);
	}

	// Confirm that carbohydrate ring conformers are loaded correctly from the database.
	void test_conformers_from_database_file_for_ring_size()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical::carbohydrates;

		TS_TRACE("Testing read_conformers_from_database_file_for_ring_size() method.");

		vector1<RingConformer> conformers =
				read_conformers_from_database_file_for_ring_size(
						"core/chemical/carbohydrates/dummy_conformers.data", 8);

		TS_ASSERT_EQUALS(conformers.size(), 3);
		TS_ASSERT_EQUALS(conformers[1].specific_name, "1F2");
		TS_ASSERT_EQUALS(conformers[2].general_name, "bar");
		TS_ASSERT_EQUALS(conformers[3].nu_angles.size(), 6);  // 2 fewer angles than the ring size should be read.
	}
};  // class CarbohydrateDatabaseIOTests
