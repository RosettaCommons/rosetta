// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/database_io.cxxtest.hh
/// @brief   Test suite for carbohydrate database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/carbohydrates/database_io.hh>

// Package header
#include <core/chemical/carbohydrates/SugarModificationsNomenclatureTable.hh>
#include <core/chemical/carbohydrates/LinkageConformers.hh>

// Project header
#include <core/types.hh>
#include <utility/excn/EXCN_Base.hh>
#include <basic/Tracer.hh>


// C++ header
#include <map>

static basic::Tracer TR("core.chemical.carbohydrates.database_io.cxxtest");

class CarbohydrateDatabaseIOTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that carbohydrate 3-letter codes and roots are loaded correctly from the database.
	void test_read_codes_and_roots_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TR << "Testing read_codes_and_roots_from_database_file() method."  << std::endl;

		map< string, RootData > map(
			read_codes_and_roots_from_database_file( "core/chemical/carbohydrates/codes_to_roots.map" ) );

		TS_ASSERT_EQUALS( map.size(), 3 );
	}

	// Confirm that carbohydrate ring sizes and their 1-letter affixes and morphemes are loaded correctly from the
	// database.
	void test_read_ring_sizes_and_morphemes_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TR << "Testing read_ring_sizes_and_morphemes_from_database_file() method."  << std::endl;

		map< core::Size, pair< char, string > > map( read_ring_sizes_and_morphemes_from_database_file(
			"core/chemical/carbohydrates/ring_size_to_morphemes.map" ) );

		TS_ASSERT_EQUALS( map.size(), 2 );
		TS_ASSERT_EQUALS( map[ 3 ].first, '\0' );  // Make sure 'X' was properly converted to a null char.
		TS_ASSERT_EQUALS( map[ 4 ].second, "ohwow" );

		// Test for bad files.
		try {
			set_throw_on_next_assertion_failure();
			read_ring_sizes_and_morphemes_from_database_file(
				"core/chemical/carbohydrates/ring_size_to_morphemes.bad_map" );
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected_message( "ERROR: read_ring_sizes_and_morphemes_from_database_file: invalid ring size; "
				"rings cannot have less than 3 atoms!" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected_message.size() ), expected_message );
		}
	}

	// TODO: Expand this test, once I've finalized what I want/need the table to contain for data.
	// Confirm that the nomenclature data table is loaded correctly from the database.
	void test_read_nomenclature_table_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TR << "Testing read_nomenclature_table_from_database_file() method."  << std::endl;

		SugarModificationsNomenclatureTable table( read_nomenclature_table_from_database_file(
			"core/chemical/carbohydrates/nomenclature.table" ) );

		TS_ASSERT_EQUALS( table.size(), 8 );  // TEMP
	}

	// Confirm that linkage data tables are loaded correctly from the database.
	void test_read_linkage_conformers_from_database_file()
	{
		using namespace std;
		using namespace core::id;
		using namespace core::chemical::carbohydrates;

		TR << "Testing read_linkage_conformers_from_database_file() method."  << std::endl;

		LinkageConformers linkages( read_linkage_conformers_from_database_file(
			"core/chemical/carbohydrates/test_linkage_data_table.tsv" ) );

		pair< string, string> key;
		utility::vector1< LinkageConformerData > conformers;

		TS_ASSERT_EQUALS( linkages.size(), 5 );  // 6 rows in the table, but 2 are for the same pair.

		key = make_pair( "alpha-evilsugar", "->3)-evilsugar" );
		conformers = linkages[ key ];
		TS_ASSERT_EQUALS( conformers.size(), 1 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_mean( phi_dihedral ), 66.6 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_sd( phi_dihedral ), 6.6 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_mean( psi_dihedral ), -66.6 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_sd( psi_dihedral ), 6.6 );
		TS_ASSERT( ! conformers[ 1 ].has_omega() );
		TS_ASSERT_EQUALS( conformers[ 1 ].n_omega(), 0 );

		key = make_pair( "beta-evilsugar", "->6)-evilsugar" );
		conformers = linkages[ key ];
		TS_ASSERT_EQUALS( conformers.size(), 1 );
		TS_ASSERT( conformers[ 1 ].has_omega() );
		TS_ASSERT_EQUALS( conformers[ 1 ].n_omega(), 1 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_mean( omega_dihedral ), 66.6 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_sd( omega_dihedral ), 6.6 );

		key = make_pair( "alpha-goodsugar", "->2)-goodsugar" );
		conformers = linkages[ key ];
		TS_ASSERT_EQUALS( conformers.size(), 2 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_mean( phi_dihedral ), 77.7 );
		TS_ASSERT_EQUALS( conformers[ 1 ].get_torsion_sd( phi_dihedral ), 7.7 );
		TS_ASSERT_EQUALS( conformers[ 2 ].get_torsion_mean( psi_dihedral ), -10.0 );
		TS_ASSERT_EQUALS( conformers[ 2 ].get_torsion_sd( psi_dihedral ), 0.0 );
		TS_ASSERT( ! conformers[ 1 ].has_omega() );
		TS_ASSERT_EQUALS( conformers[ 2 ].n_omega(), 0 );
		TS_ASSERT_EQUALS( conformers[ 1 ].population + conformers[ 2 ].population, 1.0 );

		key = make_pair( "alpha-saneshortsugar", "->9)-crazylongsugar" );
		conformers = linkages[ key ];
		TS_ASSERT_EQUALS( conformers.size(), 1 );


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
};  // class CarbohydrateDatabaseIOTests
