// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/stereochem.cxxtest.hh
/// @brief   Test suite for ensuring that the stereochemistry and default ring conformations are correct for most
/// .params files of monosaccharide residues in the database.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/rings/AxEqDesignation.hh>
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/util.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/io/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ header
#include <sstream>


static basic::Tracer TR( "core.chemical.carbohydrates.stereochem.cxxtest" );


// Type definitions
typedef utility::vector1< core::chemical::rings::AxEqDesignation > Stereochemistries;
typedef std::pair< std::string, Stereochemistries > TestKeyEntry;
typedef std::map< std::string, TestKeyEntry > TestKey;
typedef core::uint AtomIndex;


class CarbohydrateStereochemTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options( "-include_sugars" );
	}

	// Destruction
	void tearDown()
	{}

	// Return either AXIAL or EQUATORIAL from 'A' or 'E', respectively.
	core::chemical::rings::AxEqDesignation a_or_e_to_axial_or_equatorial( char letter )
	{
		using namespace core::chemical::rings;

		if ( letter == 'A' ) {
			return AXIAL;
		} else if ( letter == 'E' ) {
			return EQUATORIAL;
		} else {
			return NEITHER;
		}
	}


	// Return a map of stereochemical data to serve as a "key" for this unit test.
	TestKey read_key_from_file( std::string const & filename )
	{
		using namespace std;
		using namespace utility;


		vector1< string > const lines( io::get_lines_from_file_data( filename ) );
		TestKey key;

		for ( string const & line : lines ) {
			istringstream line_col_by_col( line );
			string code, conformer;
			char stereochem;
			Stereochemistries stereochems;

			line_col_by_col >> code >> conformer >> stereochem;
			while ( ! line_col_by_col.fail() ) {
				stereochems.push_back( a_or_e_to_axial_or_equatorial( stereochem ) );
				line_col_by_col >> stereochem;
			}

			TestKeyEntry const & key_entry( make_pair( conformer, stereochems ) );
			key[ code ] = key_entry;
		}

		return key;
	}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that all carbohydrates in the database have the correct coordinates for their virtual atoms.
	void test_virtual_atom_placement()
	{
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace chemical;

		string const filename(
			basic::database::full_name( "chemical/residue_type_sets/" + FA_STANDARD + "/" ) + "residue_types.txt" );
		vector1< string > const lines( io::get_lines_from_file_data( filename ) );
		for ( string line : lines ) {
			if ( line.substr( 0, 27 ) == "residue_types/carbohydrates"  ) {
				core::uint const start_of_name( line.find( "/to" ) );
				string const & linkage( line.substr( start_of_name + 3, 1 ) );

				string const & code(
					line.substr( start_of_name + 5, line.find( ".params" ) - ( start_of_name + 5 ) ) );
				if ( code == "alpha-Daup" ) { continue; }  // TEMP
				TR << "Testing: ->" << linkage << ")-" << code << endl;
				pose::Pose pose;
				make_pose_from_saccharide_sequence( pose, "->" + linkage + ")-" + code );
				//pose.dump_pdb( "/home/labonte/Dropbox/Transfer/Sugars/" + linkage + "-" + code + ".pdb" );  // TEMP
				conformation::Residue const & res( pose.residue( 1 ) );
				TS_ASSERT( res.is_carbohydrate() );
				carbohydrates::CarbohydrateInfoCOP info( res.carbohydrate_info() );

				TS_ASSERT_DELTA( res.xyz( info->anomeric_carbon_index() ).x(),
					res.xyz( info->virtual_anomeric_carbon_index() ).x(), 0.05 );
				TS_ASSERT_DELTA( res.xyz( info->anomeric_carbon_index() ).y(),
					res.xyz( info->virtual_anomeric_carbon_index() ).y(), 0.05 );
				TS_ASSERT_DELTA( res.xyz( info->anomeric_carbon_index() ).z(),
					res.xyz( info->virtual_anomeric_carbon_index() ).z(), 0.05 );

				TS_ASSERT_DELTA( res.xyz( info->cyclic_oxygen_index() ).x(),
					res.xyz( info->virtual_cyclic_oxygen_index() ).x(), 0.05 );
				TS_ASSERT_DELTA( res.xyz( info->cyclic_oxygen_index() ).y(),
					res.xyz( info->virtual_cyclic_oxygen_index() ).y(), 0.05 );
				TS_ASSERT_DELTA( res.xyz( info->cyclic_oxygen_index() ).z(),
					res.xyz( info->virtual_cyclic_oxygen_index() ).z(), 0.05 );
			}
		}
	}

	// Confirm that carbohydrates in the database have the correct stereochemistry and default ring conformers.
	void test_rings_and_stereochemistries()
	{
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace chemical;
		using namespace rings;
		using namespace conformation;

		TestKey const key( read_key_from_file( "core/chemical/carbohydrates/stereochem_test.key" ) );
		TS_ASSERT_EQUALS( key.size(), 88 );  // 16 aldopentoses & 32 aldohexoses & 40 common modified sugars

		pose::Pose pose;
		GlobalResidueTypeSetCOP type_set( new GlobalResidueTypeSet(
			FA_STANDARD, basic::database::full_name( "chemical/residue_type_sets/" + FA_STANDARD + "/" ) ) );

		for ( auto const & code_key_entry_pair : key ) {
			string const & code( code_key_entry_pair.first );
			TestKeyEntry const & key_entry( code_key_entry_pair.second );
			string const & conformer( key_entry.first );
			Stereochemistries const & stereochemistries( key_entry.second );

			for ( core::uint linkage( 1 );
					linkage <= chemical::carbohydrates::CarbohydrateInfo::MAX_C_SIZE_LIMIT; ++linkage ) {
				string const name( "->" + to_string( linkage ) + ")-" + code );

				// Check that the residue exists in the database.
				if ( type_set->has_name( name ) ) {
					TR << "Testing: " << name << endl;

					// Make a disaccharide, to test both reducing and non-reducing ends.
					// TODO: To test ketoses, this line needs to be changed.
					make_pose_from_saccharide_sequence( pose, name + "-(1" + name );
					Residue const & reducing_end( pose.residue( 1 ) );
					Residue const & non_reducing_end( pose.residue( 2 ) );

					// Check that the ring conformers are as expected.
					TR << "  Testing ring with nu angles: ";
					for ( Angle nu : reducing_end.nus() ) {
						TR << nu << " ";
					}
					TR << endl;
					TS_ASSERT_EQUALS( reducing_end.ring_conformer( 1 ).specific_name, conformer );
					TS_ASSERT_EQUALS( non_reducing_end.ring_conformer( 1 ).specific_name, conformer );

					// Check that the stereochemistries are as expected.
					AtomIndices const & ring_atoms( reducing_end.type().ring_atoms( 1 ) );
					for ( core::uint i( 1 ); i <= stereochemistries.size(); ++i ) {
						AtomIndex const ring_atom( ring_atoms[ i ] );

						TR << "  Testing stereochemistry at ";
						TR << reducing_end.atom_name( ring_atom ) << ' ';

						AtomIndices const substituents( reducing_end.get_substituents_to_ring_atom( ring_atom ) );

						// Hydrogen indices must be stored separately, because the indices change between internal and
						// terminal residues, since heavy atoms are added in the latter case.
						AtomIndices const reducing_end_hydrogens(
							reducing_end.get_hydrogens_bonded_to_ring_atom( ring_atom ) );
						AtomIndices const non_reducing_end_hydrogens(
							non_reducing_end.get_hydrogens_bonded_to_ring_atom( ring_atom ) );

						if ( substituents.size() == 0 ) {
							TR << endl;
							TS_ASSERT_EQUALS( reducing_end_hydrogens.size(), 2 );
							TS_ASSERT_EQUALS( non_reducing_end_hydrogens.size(), 2 );

							// If there is no substituent, this is a deoxy sugar.
							// In that case, we just confirm that the designation is N/A and move on.
							TS_ASSERT_EQUALS( stereochemistries.at( i ), NEITHER );
							continue;
						}

						// Otherwise, we have at least 1 substituent.
						AtomIndex const substituent( substituents[ 1 ] );

						// TODO: I need to check if there are any whacky dehydro sugars.
						AtomIndex const reducing_end_hydrogen( reducing_end_hydrogens[ 1 ] );
						AtomIndex const non_reducing_end_hydrogen( non_reducing_end_hydrogens[ 1 ] );

						TR << "using " << reducing_end.atom_name( substituent ) << endl;
						TS_ASSERT_EQUALS(
							is_atom_axial_or_equatorial( reducing_end, substituent ), stereochemistries.at( i ) );

						// If this happens to be the anomeric carbon, we need to test the UPPER connect of a parent
						// residue to confirm that the LOWER connection was placed properly.
						// Since the atoms are on two distinct residues, we need to use coordinates.
						if ( ring_atom == non_reducing_end.carbohydrate_info()->anomeric_carbon_index() ) {
							Coords const linkage_atom_coords( reducing_end.xyz( reducing_end.upper_connect_atom() ) );
							Coords const anomeric_carbon_coords( non_reducing_end.xyz( ring_atom ) );
							Size const n_ring_atoms( ring_atoms.size() );
							vector1< Coords > ring_atom_coords( n_ring_atoms );
							for ( core::uint j( 1 ); j <= n_ring_atoms; ++j ) {
								ring_atom_coords[ j ] = non_reducing_end.xyz( ring_atoms[ j ] );
							}

							TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
								linkage_atom_coords, anomeric_carbon_coords, ring_atom_coords ),
								stereochemistries.at( i ) );
						}

						TR << "  Testing stereochemistry of the hydrogens." << endl;
						TS_ASSERT_EQUALS( is_atom_axial_or_equatorial( reducing_end, reducing_end_hydrogen ),
							opposite_designation( stereochemistries.at( i ) ) );
						TS_ASSERT_EQUALS( is_atom_axial_or_equatorial( non_reducing_end, non_reducing_end_hydrogen ),
							opposite_designation( stereochemistries.at( i ) ) );

						TR << "  Testing placement of the hydrogens." << endl;
						for ( AtomIndex ring_atom : ring_atoms ) {
							TS_ASSERT( reducing_end.xyz( reducing_end_hydrogen ).distance( reducing_end.xyz( ring_atom ) )
								> 1.0 );  // If it's closer than 1 angstrom, it's certainly incorrectly placed.
							TS_ASSERT( non_reducing_end.xyz( non_reducing_end_hydrogen ).distance( non_reducing_end.xyz( ring_atom ) )
								> 1.0 );  // If it's closer than 1 angstrom, it's certainly incorrectly placed.
						}
					}
				} else {
					TR.Debug << "Residue " << name << " not found in database; skipping." << endl;
				}
			}
		}
	}

};  // class CarbohydrateStereochemTests
