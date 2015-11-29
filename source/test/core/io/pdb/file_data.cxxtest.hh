// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/io/pdb/file_data.cxxtest.hh
/// @brief   Test suite for FileData methods and utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/file_data.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>


class FileDataTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Input Methods //////////////////////////////////////////////////////////
	// Confirm that LINK records are stored properly.
	void test_store_link_record()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io::pdb;
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1");
		TS_TRACE( "Testing store_link_record() method of FileData." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
		string const sample_pdb_lines(
			"LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98  \n"
			"LINK         OG1 THR A  51                NA    NA A6001     1555   1555  2.72  \n"
			"LINK         OD2 ASP A  66                NA    NA A6001     1555   1555  2.72  \n"
			"LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93  \n"
			"LINK         C21 2EG A   7                 C22 2EG B  19     1555   1555  1.56  \n" );

		utility::vector1< Record > records( PDB_DReader::parse( sample_pdb_lines ) );

		FileData fd;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 5 );
		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			fd.store_link_record( records[ i ] );
		}

		map< string, vector1< LinkInformation > > link_map( fd.link_map );
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

	// Confirm that SSBOND records are stored properly.
	void test_store_ssbond_record()
	{
		using namespace std;
		using namespace utility;
		using namespace core::io::pdb;
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1");
		TS_TRACE( "Testing store_ssbond_record() method of FileData." );

		// These example lines are taken directly from the given examples at
		// http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
		string const sample_pdb_lines(
			"SSBOND   1 CYS A    6    CYS A  127                          1555   1555  2.03  \n"
			"SSBOND   2 CYS A   30    CYS A  115                          1555   1555  2.07  \n"
			"SSBOND   3 CYS A   64    CYS A   80                          1555   1555  2.06  \n"
			"SSBOND   4 CYS A   76    CYS A   94                          1555   1555  2.04  \n" );

		utility::vector1< Record > records( PDB_DReader::parse( sample_pdb_lines ) );

		FileData fd;

		core::Size const n_records( records.size() );
		TS_ASSERT_EQUALS( n_records, 4 );
		for ( core::uint i( 1 ); i <= n_records; ++i ) {
			fd.store_ssbond_record( records[ i ] );
		}

		map< string, vector1< SSBondInformation > > ssbond_map( fd.ssbond_map );
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

	// Confirm that HETNAM records are stored properly.
	void test_store_heterogen_names()
	{
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1");
		TS_ASSERT( true );
	}

	// Output Methods /////////////////////////////////////////////////////////
	// Confirm that LinkInformation and SSBondInformation data are created properly from a Pose.
	void test_get_connectivity_annotation_info()
	{
		using namespace core::io::pdb;
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1");

		TS_TRACE( "Testing get_connectivity_annotation_info() method of FileData." );

		core::pose::Pose pose;

		// 1BH4 is circulin A, a cyclic peptide with three disulfides.
		//core::import_pose::pose_from_pdb( pose, "core/io/pdb/1BH4.pdb" );

		// 4TTL is a cyclic peptide with two disulfides.
		core::import_pose::pose_from_pdb( pose, "core/io/pdb/4TTL.pdb" );

		//core::import_pose::pose_from_pdb( pose, "core/io/pdb/1e68_link.pdb" );

		FileData fd;

		TS_ASSERT_EQUALS( fd.link_map.size(),   0 );
		fd.get_connectivity_annotation_info( pose );

		TS_ASSERT_EQUALS( fd.ssbond_map.size(), 2 );
		TS_ASSERT_EQUALS( fd.link_map.size(),   1 );

		// The following lines came directly from 1BH4:
		//SSBOND   1 CYS A    1    CYS A   17                          1555   1555  2.02
		//SSBOND   2 CYS A    5    CYS A   19                          1555   1555  2.02
		//SSBOND   3 CYS A   10    CYS A   24                          1555   1555  2.02
		//LINK         N   CYS A   1                 C   PRO A  30     1555   1555  1.31

		// The following lines came directly from 4TTL:
		//SSBOND   1 CYS A    2    CYS A    8                          1555   1555  2.05
		//SSBOND   2 CYS A    3    CYS A   16                          1555   1555  2.03
		//LINK         N   GLY A   1                 C   GLY A  22     1555   1555  1.34

		// The following lines come directly from 1e68_link.pdb:
		//LINK         N   MET A   1                 C   TRP A  70     1555   1555  1.33

		core::pose::Pose pose2;
		core::import_pose::pose_from_pdb( pose2, "core/io/pdb/1e68_link.pdb" );
		FileData fd2;

		TS_ASSERT_EQUALS( fd2.link_map.size(),   0 );
		fd2.get_connectivity_annotation_info( pose2 );

		TS_ASSERT_EQUALS( fd2.ssbond_map.size(), 0 );
		TS_ASSERT_EQUALS( fd2.link_map.size(),   1 );
	}
	
	/// @brief Tests PDB import when the -extra_res_fa flag is used to specify a noncanonical residue.
	///
	void test_extra_res_fa_flag() {
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1 -extra_res_fa core/io/pdb/test.params");
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "core/io/pdb/extra_res_pose.pdb" );
		TS_ASSERT_EQUALS( pose.n_residue(), 181 );
	}

};  // class FileDataTests

// Sandbox ////////////////////////////////////////////////////////////////////


