// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/io/pdb/file_data.cxxtest.hh
/// @brief   Test suite for StructFileRep methods and utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>


class StructFileRepTests : public CxxTest::TestSuite {
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


	/// @brief Tests PDB import when the -extra_res_fa flag is used to specify a noncanonical residue.
	void test_extra_res_fa_flag() {
		// This has to be the first one called because we don't re-initialize core and we really need to get this extra res
		// in the RTS...
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1 -extra_res_fa core/io/pdb/test.params");
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/io/pdb/extra_res_pose.pdb" , core::import_pose::PDB_file);
		TS_ASSERT_EQUALS( pose.n_residue(), 181 );
	}


	// Output Methods /////////////////////////////////////////////////////////
	// Confirm that LinkInformation and SSBondInformation data are created properly from a Pose.
	void test_get_connectivity_annotation_info()
	{
		using namespace core::io;
		using namespace core::io::pdb;
		core_init_with_additional_options( "-obey_ENDMDL -read_pdb_link_records -constraints_from_link_records -cst_weight 1");

		TS_TRACE( "Testing get_connectivity_annotation_info() method of StructFileRep." );

		core::pose::Pose pose;

		// 1BH4 is circulin A, a cyclic peptide with three disulfides.
		//core::import_pose::pose_from_file( pose, "core/io/pdb/1BH4.pdb" , core::import_pose::PDB_file);

		// 4TTL is a cyclic peptide with two disulfides.
		core::import_pose::pose_from_file( pose, "core/io/pdb/4TTL.pdb" , core::import_pose::PDB_file);

		//core::import_pose::pose_from_file( pose, "core/io/pdb/1e68_link.pdb" , core::import_pose::PDB_file);

		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->link_map().size(),   0 );
		pose_to_sfr.get_connectivity_annotation_info( pose );

		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->ssbond_map().size(), 2 );
		TS_ASSERT_EQUALS( pose_to_sfr.sfr()->link_map().size(),   1 );

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
		core::import_pose::pose_from_file( pose2, "core/io/pdb/1e68_link.pdb" , core::import_pose::PDB_file);
		core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr2;

		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->link_map().size(),   0 );
		pose_to_sfr2.get_connectivity_annotation_info( pose2 );

		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->ssbond_map().size(), 0 );
		TS_ASSERT_EQUALS( pose_to_sfr2.sfr()->link_map().size(),   1 );
	}

};  // class StructFileRepTests

// Sandbox ////////////////////////////////////////////////////////////////////


