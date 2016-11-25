// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/carbohydrates/CarbohydrateInfoManager.cxxtest.hh
/// @brief   Test suite for CarbohydrateInfoManager
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/LinkageConformers.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("core.chemical.carbohydrates.CarbohydrateInfoManager.cxxtest");

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
		using namespace core::id;

		TR << "Testing access to linkage conformer data..."  << std::endl;

		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-evilsugar", "alpha-evilsugar" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->3)-alpha-evilsugar-" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->X)-alpha-evilsugar-" ) );
		TS_ASSERT( CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-beta-evilsugar-", "->3)-alpha-evilsugar-" ) );
		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics( "->3)-alpha-evilsugar-", "->3)-beta-evilsugar-" ) );
		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics( "->X)-alpha-evilsugar-", "->3)-alpha-evilsugar-" ) );
		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics( "->X)-alpha-evilsugar-", "->3)-alpha-evilsugar-" ) );

		TS_ASSERT( ! CarbohydrateInfoManager::pair_has_linkage_statistics(  "alpha-saneshortsugar", "->9)-crazylongsugar") );
		TS_ASSERT_EQUALS(
			CarbohydrateInfoManager::linkages_from_pair( "->2)-goodsugar", "alpha-goodsugar" ).size(), 2 );



		utility::vector1< LinkageConformerData > const & conformers = CarbohydrateInfoManager::linkages_from_pair( "->9)-crazylongsugar", "alpha-saneshortsugar");

		LinkageConformerData data = conformers[ 1 ];

		//Test Phi and Psi
		//12.0 4.0 67.0 9.0 12.0 4.0 65.0 9.0 66.0 4.0 67.0 9.0
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( phi_dihedral)), 12);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( phi_dihedral)), 4);

		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( psi_dihedral)), 67);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( psi_dihedral)), 9);

		//Test Omegas
		TS_ASSERT_EQUALS(
			conformers.size(), 1 );

		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( omega_dihedral, 1)), 12);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( omega_dihedral, 1)), 4);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( omega_dihedral, 2)), 65);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( omega_dihedral, 2)), 9);

		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( omega_dihedral, 3)), 66);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( omega_dihedral, 3)), 5);

		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_mean( omega_dihedral, 4)), 67);
		TS_ASSERT_EQUALS(
			core::Size(data.get_torsion_sd( omega_dihedral, 4)), 10);


	}

};  // class CarbohydrateInfoManagerTests
