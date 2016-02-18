// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/chemical/carbohydrates/CarbohydrateInfoManager.cxxtest.hh
/// @brief   Test suite for CarbohydrateInfoManager
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

// Utility headers
#include <utility/vector1.hh>


class CarbohydrateInfoManagerTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options(
				"-linkage_conformer_data_file core/chemical/carbohydrates/linkage_data_table.tsv" );
	}
	
	// Destruction
	void tearDown()
	{}
	

public:  // Tests /////////////////////////////////////////////////////////////
	void test_linkage_conformer_data_access()
	{
		using namespace core::chemical::carbohydrates;

		TS_TRACE( "Testing access to linkage conformer data..." );

		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-evilsugar", "alpha-evilsugar" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->3)-alpha-evilsugar-" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->X)-alpha-evilsugar-" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-beta-evilsugar-", "->3)-alpha-evilsugar-" ) );
		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->3)-beta-evilsugar-" ) );
		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics( "->X)-alpha-evilsugar-", "->3)-alpha-evilsugar-" ) );
		
		TS_ASSERT_EQUALS(
				CarbohydrateInfoManager::linkages_from_pair( "->2)-goodsugar", "alpha-goodsugar" ).size(), 2 );
	}

};  // class CarbohydrateInfoManagerTests
