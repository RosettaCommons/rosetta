// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/io/pdb/record_def_io.cxxtest.hh
/// @brief   Test suite for database loading of the reference table of PDB record definitions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/pdb_reader.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ header
#include <map>


class PDBReaderTests : public CxxTest::TestSuite {
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
	// Confirm that SEQRES records are stored properly.
	void test_store_chain_sequence_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_chain_sequence_record_in_sfr() method." );

		// These example lines are taken almost directly from the given examples at
		// http://http://www.wwpdb.org/documentation/file-format-content/format33/sect3.html
		string const sample_pdb_lines(
			"SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU          \n"
			"SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN                              \n"
			"SEQRES   1 B   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU          \n"
			"SEQRES   2 B   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR          \n"
			"SEQRES   3 B   30  THR PRO LYS ALA                                              \n"
			"SEQRES   1 C   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU          \n"
			"SEQRES   2 C   21  TYR GLN LEU GLU ASN TYR CYS ASN                              \n"
			"SEQRES   1 D   26  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU          \n"
			"SEQRES   2 D   26  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR          \n"
			"SEQRES   1 E    8   DA  DA  DC  DC  DG  DG  DT  DT                              \n"
			"SEQRES   1 F    8   DA  DA  DC  DC  DG  DG  DT  DT                              \n"
			"SEQRES   1 X   39    U   C   C   C   C   C   G   U   G   C   C   C   A          \n"
			"SEQRES   2 X   39    U   A   G   C   G   G   C   G   U   G   G   A   A          \n"
			"SEQRES   3 X   39    C   C   A   C   C   C   G   U   U   C   C   C   A          \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 14 );

		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_chain_sequence_record_in_sfr( records[ i ], sfr );
		}

		map< char, vector1< string > > sequences_map( sfr.chain_sequences() );

		TS_ASSERT_EQUALS( sequences_map.size(), 7 );  // 7 different chains

		TS_ASSERT( sequences_map.count( 'A' ) );
		TS_ASSERT( sequences_map.count( 'B' ) );
		TS_ASSERT( sequences_map.count( 'C' ) );
		TS_ASSERT( sequences_map.count( 'D' ) );
		TS_ASSERT( sequences_map.count( 'E' ) );
		TS_ASSERT( sequences_map.count( 'F' ) );
		TS_ASSERT( sequences_map.count( 'X' ) );

		TS_ASSERT_EQUALS( sequences_map[ 'A' ].size(), 21 );
		TS_ASSERT_EQUALS( sequences_map[ 'B' ].size(), 30 );
		TS_ASSERT_EQUALS( sequences_map[ 'C' ].size(), 21 );
		TS_ASSERT_EQUALS( sequences_map[ 'D' ].size(), 26 );
		TS_ASSERT_EQUALS( sequences_map[ 'E' ].size(), 8 );
		TS_ASSERT_EQUALS( sequences_map[ 'F' ].size(), 8 );
		TS_ASSERT_EQUALS( sequences_map[ 'X' ].size(), 39 );

		string sequenceA, sequenceB, sequenceE, sequenceX;
		for ( core::uint i( 1 ); i <= 21; ++i ) {
			sequenceA += ( sequences_map[ 'A' ][ i ] + " " );
		}
		TS_ASSERT_EQUALS( sequenceA, "GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU "
			"TYR GLN LEU GLU ASN TYR CYS ASN " );
		for ( core::uint i( 1 ); i <= 30; ++i ) {
			sequenceB += ( sequences_map[ 'B' ][ i ] + " " );
		}
		TS_ASSERT_EQUALS( sequenceB, "PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU "
			"ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR "
			"THR PRO LYS ALA ");
		for ( core::uint i( 1 ); i <= 8; ++i ) {
			sequenceE += ( sequences_map[ 'E' ][ i ] + " " );
		}
		TS_ASSERT_EQUALS( sequenceE, " DA  DA  DC  DC  DG  DG  DT  DT " );
		for ( core::uint i( 1 ); i <= 39; ++i ) {
			sequenceX += ( sequences_map[ 'X' ][ i ] + " " );
		}
		TS_ASSERT_EQUALS( sequenceX, "  U   C   C   C   C   C   G   U   G   C   C   C   A "
			"  U   A   G   C   G   G   C   G   U   G   G   A   A "
			"  C   C   A   C   C   C   G   U   U   C   C   C   A " );
	}

	// Confirm that MODRES records are stored properly.
	void test_store_mod_res_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_mod_res_record_in_sfr() method." );

		// These example lines are taken almost directly from the given examples at
		// http://http://www.wwpdb.org/documentation/file-format-content/format33/sect3.html
		string const sample_pdb_lines(
			"MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE                                    \n"
			"MODRES 1IL2 1MG D 1937    G  1N-METHYLGUANOSINE-5'-MONOPHOSPHATE                   \n"
			"MODRES 4ABC MSE B   32  MET  SELENOMETHIONINE                                      \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 3 );

		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_mod_res_record_in_sfr( records[ i ], sfr );
		}

		map< string, ModifiedResidueInformation > modres_map( sfr.modres_map() );
		TS_ASSERT_EQUALS( modres_map.size(), 3 );

		TS_ASSERT( modres_map.count( "  74 A" ) );
		TS_ASSERT( modres_map.count( "1937 D" ) );
		TS_ASSERT( modres_map.count( "  32 B" ) );

		TS_ASSERT_EQUALS( modres_map[ "  74 A" ].comment, "GLYCOSYLATION SITE" );
		TS_ASSERT_EQUALS( modres_map[ "1937 D" ].stdRes, "  G" );
		TS_ASSERT_EQUALS( modres_map[ "  32 B" ].iCode, ' ' );
		TS_ASSERT_EQUALS( modres_map[ "  32 B" ].seqNum, 32 );
		TS_ASSERT_EQUALS( modres_map[ "  32 B" ].chainID, 'B' );
	}

	// Confirm that HETNAM records are stored properly.
	void test_store_heterogen_name_record_in_sfr_and_store_base_residue_type_name_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_heterogen_name_record_in_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html
		// with the addition of a single "Rosetta-format" line.
		string const sample_pdb_lines(
				"HETNAM     NAG N-ACETYL-D-GLUCOSAMINE                                          \n"
				"HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE                 \n"
				"HETNAM  2  SAD DINUCLEOTIDE                                                    \n"
				"HETNAM     UDP URIDINE-5'-DIPHOSPHATE                                          \n"
				"HETNAM     UNX UNKNOWN ATOM OR ION                                             \n"
				"HETNAM     UNL UNKNOWN LIGAND                                                  \n"
				"HETNAM     B3P 2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-                \n"
				"HETNAM   2 B3P  PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL                  \n"
				"HETNAM     Krp X  13Z Kryptonite, which will kill Superman -- Bwahaha!         \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 9 );

		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_heterogen_name_record_in_sfr( records[ i ], sfr );
		}

		map< string, string > names_map( sfr.heterogen_names() );
		TS_ASSERT_EQUALS( names_map.size(), 7 );

		TS_ASSERT( names_map.count( "NAG" ) );
		TS_ASSERT( names_map.count( "SAD" ) );
		TS_ASSERT( names_map.count( "UDP" ) );
		TS_ASSERT( names_map.count( "UNX" ) );
		TS_ASSERT( names_map.count( "UNL" ) );
		TS_ASSERT( names_map.count( "B3P" ) );
		TS_ASSERT( names_map.count( "Krp" ) );

		TS_ASSERT_EQUALS( names_map[ "NAG" ], "N-ACETYL-D-GLUCOSAMINE" );
		TS_ASSERT_EQUALS( names_map[ "SAD" ], "BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE DINUCLEOTIDE" );
		TS_ASSERT_EQUALS( names_map[ "UDP" ], "URIDINE-5'-DIPHOSPHATE" );
		TS_ASSERT_EQUALS( names_map[ "UNX" ], "UNKNOWN ATOM OR ION" );
		TS_ASSERT_EQUALS( names_map[ "UNL" ], "UNKNOWN LIGAND" );
		TS_ASSERT_EQUALS( names_map[ "B3P" ],
				"2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL" );
		TS_ASSERT_EQUALS( names_map[ "Krp" ], "X  13Z Kryptonite, which will kill Superman -- Bwahaha!" );


		TS_TRACE( "Testing (indirectly) store_base_residue_type_name_in_sfr method, "
				"which is called from within store_heterogen_name_record_in_sfr()...");

		map< string, pair< string, string > > base_names( sfr.residue_type_base_names() );

		// Each record is treated as separate for base names.
		// However, because the resID (See below.) for both "UNKNOWN ATOM OR
		// ION" and "UNKNOWN LIGAND" are the same, the latter overwrites the
		// former, giving us 8 records instead of 9.
		TS_ASSERT_EQUALS( base_names.size(), 8 );

		// Each key here is a 6-character "resID" used internally by Rosetta during Pose building/unbuilding.
		TS_ASSERT( base_names.count( "-ACETN" ) );  // a meaningless code, since this record was a PDB-format record
		TS_ASSERT( base_names.count( "  13ZX" ) );  // residue 13Z of chain X

		TS_ASSERT_EQUALS( base_names[ "-ACETN" ].first, "NAG" );
		TS_ASSERT_EQUALS( base_names[ "-ACETN" ].second, "L-D-GLUCOSAMINE" );  // meaningless junk

		TS_ASSERT_EQUALS( base_names[ "  13ZX" ].first, "Krp" );
		TS_ASSERT_EQUALS( base_names[ "  13ZX" ].second, "Kryptonite" );  // Any word after a comma is ignored.
	}

	// Confirm that HETSYM records are stored properly.
	void test_store_heterogen_synonym_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_heterogen_synonym_record_in_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html
		string const sample_pdb_lines(
				"HETSYN     HV5 3-METHYL-L-VALINE                                               \n"
				"HETSYN     AB1 ABT-378; LOPINAVIR                                              \n"
				"HETSYN     CMP CYCLIC AMP; CAMP                                                \n"
				"HETSYN     TRS TRIS  BUFFER;                                                   \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 4 );

		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_heterogen_synonym_record_in_sfr( records[ i ], sfr );
		}

		map< string, vector1< string > > synonyms_map( sfr.heterogen_synonyms() );
		TS_ASSERT_EQUALS( synonyms_map.size(), 4 );

		TS_ASSERT( synonyms_map.count( "HV5" ) );
		TS_ASSERT( synonyms_map.count( "AB1" ) );
		TS_ASSERT( synonyms_map.count( "CMP" ) );
		TS_ASSERT( synonyms_map.count( "TRS" ) );

		TS_ASSERT_EQUALS( synonyms_map[ "HV5" ][ 1 ], "3-METHYL-L-VALINE" );
		TS_ASSERT_EQUALS( synonyms_map[ "AB1" ][ 2 ], "LOPINAVIR" );
		TS_ASSERT_EQUALS( synonyms_map[ "CMP" ][ 1 ], "CYCLIC AMP" );
		TS_ASSERT_EQUALS( synonyms_map[ "TRS" ][ 1 ], "TRIS  BUFFER" );
	}

	// Confirm that FORMUL records are stored properly.
	void test_store_formula_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_formula_record_in_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html
		string const sample_pdb_lines(
			"FORMUL   3   MG    2(MG 2+)                                                    \n"
			"FORMUL   5  SO4    6(O4 S 2-)                                                  \n"
			"FORMUL  13  HOH   *360(H2 O)                                                   \n"
			"FORMUL   3  NAP    2(C21 H28 N7 O17 P3)                                        \n"
			"FORMUL   4  FOL    2(C19 H19 N7 O6)                                            \n"
			"FORMUL   5  1PE    C10 H22 O6                                                  \n"
			"FORMUL   2  NX5    C14 H10 O2 CL2 S                                            \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 7 );

		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_formula_record_in_sfr( records[ i ], sfr );
		}

		map< string, string > formulae( sfr.heterogen_formulae() );
		TS_ASSERT_EQUALS( formulae.size(), 7 );

		TS_ASSERT( formulae.count( " MG" ) );
		TS_ASSERT( formulae.count( "SO4" ) );
		TS_ASSERT( formulae.count( "HOH" ) );
		TS_ASSERT( formulae.count( "NAP" ) );
		TS_ASSERT( formulae.count( "FOL" ) );
		TS_ASSERT( formulae.count( "1PE" ) );
		TS_ASSERT( formulae.count( "NX5" ) );

		TS_ASSERT_EQUALS( formulae[ "HOH" ], "*360(H2 O)" );
		TS_ASSERT_EQUALS( formulae[ "NAP" ], " 2(C21 H28 N7 O17 P3)" );
	}

	// Confirm that SSBOND records are stored properly.
	void test_store_ssbond_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_ssbond_record_in_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
		string const sample_pdb_lines(
			"SSBOND   1 CYS A    6    CYS A  127                          1555   1555  2.03  \n"
			"SSBOND   2 CYS A   30    CYS A  115                          1555   1555  2.07  \n"
			"SSBOND   3 CYS A   64    CYS A   80                          1555   1555  2.06  \n"
			"SSBOND   4 CYS A   76    CYS A   94                          1555   1555  2.04  \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 4 );
		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_ssbond_record_in_sfr( records[ i ], sfr );
		}

		map< string, vector1< SSBondInformation > > ssbond_map( sfr.ssbond_map() );
		TS_ASSERT_EQUALS( ssbond_map.size(), 4 );

		TS_ASSERT( ssbond_map.count( "  30 A" ) );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ].size(), 1 );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].resName1, "CYS" );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].chainID1, 'A' );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].resSeq1, 30 );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].iCode1, ' ' );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].resName2, "CYS" );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].chainID2, 'A' );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].resSeq2, 115 );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].iCode2, ' ' );
		TS_ASSERT_EQUALS( ssbond_map[ "  30 A" ][ 1 ].length, 2.07 );

		TS_ASSERT( ssbond_map.count( "  64 A" ) );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ].size(), 1 );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].resName1, "CYS" );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].chainID1, 'A' );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].resSeq1, 64 );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].iCode1, ' ' );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].resName2, "CYS" );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].chainID2, 'A' );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].resSeq2, 80 );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].iCode2, ' ' );
		TS_ASSERT_EQUALS( ssbond_map[ "  64 A" ][ 1 ].length, 2.06 );
	}

	// Confirm that LINK records are stored properly.
	void test_store_link_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_link_recordin_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
		string const sample_pdb_lines(
			"LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98  \n"
			"LINK         OG1 THR A  51                NA    NA A6001     1555   1555  2.72  \n"
			"LINK         OD2 ASP A  66                NA    NA A6001     1555   1555  2.72  \n"
			"LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93  \n"
			"LINK         C21 2EG A   7                 C22 2EG B  19     1555   1555  1.56  \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 5 );
		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_link_record_in_sfr( records[ i ], sfr );
		}

		map< string, vector1< LinkInformation > > link_map( sfr.link_map() );
		TS_ASSERT_EQUALS( link_map.size(), 5 );

		TS_ASSERT( link_map.count( "  49 A" ) );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ].size(), 1 );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].name1, " O  " );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].resName1, "GLY" );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].chainID1, 'A' );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].resSeq1, 49 );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].iCode1, ' ' );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].name2, "NA  " );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].resName2, " NA" );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].chainID2, 'A' );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].resSeq2, 6001 );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].iCode2, ' ' );
		TS_ASSERT_EQUALS( link_map[ "  49 A" ][ 1 ].length, 2.98 );

		TS_ASSERT( link_map.count( "   7 A" ) );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ].size(), 1 );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].name1, " C21" );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].resName1, "2EG" );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].chainID1, 'A' );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].resSeq1, 7 );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].iCode1, ' ' );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].name2, " C22" );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].resName2, "2EG" );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].chainID2, 'B' );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].resSeq2, 19 );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].iCode2, ' ' );
		TS_ASSERT_EQUALS( link_map[ "   7 A" ][ 1 ].length, 1.56 );
	}

	// Confirm that CISPEP records are stored properly.
	void test_store_cis_peptide_record_in_sfr()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_cis_peptide_record_in_sfr() method." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
		string const sample_pdb_lines(
			"CISPEP   1 SER A   58    GLY A   59          0        20.91                    \n"
			"CISPEP   1 GLY A  116    GLY A  117          0        18.50                    \n"
			"CISPEP   1 MET A    1    SER A    2          0        -3.69                    \n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		StructFileRep sfr;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 3 );
		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			store_cis_peptide_record_in_sfr( records[ i ], sfr );
		}

		map< string, CisPeptideInformation > cispep_map( sfr.cispep_map() );
		TS_ASSERT_EQUALS( cispep_map.size(), 3 );

		TS_ASSERT( cispep_map.count( "  58 A" ) );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].pep1, "SER" );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].chainID1, 'A' );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].seqNum1, 58 );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].icode1, ' ' );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].pep2, "GLY" );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].chainID2, 'A' );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].seqNum2, 59 );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].icode2, ' ' );
		TS_ASSERT_EQUALS( cispep_map[ "  58 A" ].measure, 20.91 );

		TS_ASSERT( cispep_map.count( " 116 A" ) );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].pep1, "GLY" );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].chainID1, 'A' );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].seqNum1, 116 );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].icode1, ' ' );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].pep2, "GLY" );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].chainID2, 'A' );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].seqNum2, 117 );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].icode2, ' ' );
		TS_ASSERT_EQUALS( cispep_map[ " 116 A" ].measure, 18.5 );
	}

	// Confirm that UNKNOW records are stored properly.
	void test_store_unknown_records_in_sfr() {
		using namespace std;
		using namespace utility;
		using namespace core::io;
		using namespace core::io::pdb;

		TS_TRACE( "Testing store_unknown_records_in_sfr() method." );

		string const sample_pdb_lines(
			"ATOM   8914  OXT GLY G 472     -31.771  31.214 -89.955  1.00  0.00           O  \n"
			"ENDMDL                                                                          \n"
			"##Begin comments##\n"
			"ACCEPT LOG 1 -1 -418.883 NA\n"
			"GRAFT LOG 1 DATA GRAFT_CLOSURE 1 input_pdbs/pareto_4JAN_CH103_GP120_renum_0001.pdb L3 L3-9-2 4fzeL_L3_0001\n"
			"REMARK H1 CLUSTER H1-13-4 29.3958\n"
			"##End comments##\n"
			"label fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama omega fa_dun p_aa_pp ref total\n"
			"weights 0.8 0.44 0.75 0.004 0.7 1 1.17 1.17 1.17 1.1 1 0.2 0.5 0.56 0.32 1 NA\n"
			"pose -2181.57 216.307 1195.77 5.07317 -229.036 101.933 -52.4809 -228.717 -76.359 -78.7517 -3.85772 -19.3249 74.1669 471.109 -89.7754 19.6818 -875.835\n" );

		vector1< Record > records( create_records_from_pdb_file_contents( sample_pdb_lines ) );

		TS_ASSERT_EQUALS( records.size(), 10 );

		StructFileReaderOptions options;
		options.set_pdb_comments( true );
		StructFileRep sfr( create_sfr_from_pdb_records( records, options ) );

		map< string, string > comments( sfr.pdb_comments() );

		// REMARK is a .pdb standard record and should have been loaded earlier.
		TS_ASSERT_EQUALS( comments.size(), 2 );

		TS_ASSERT_EQUALS( comments[ "ACCEPT" ], "LOG 1 -1 -418.883 NA" );
		TS_ASSERT_EQUALS( comments[ "GRAFT" ],
			"LOG 1 DATA GRAFT_CLOSURE 1 input_pdbs/pareto_4JAN_CH103_GP120_renum_0001.pdb L3 L3-9-2 4fzeL_L3_0001" );
	}
};  // class PDBReaderTests
