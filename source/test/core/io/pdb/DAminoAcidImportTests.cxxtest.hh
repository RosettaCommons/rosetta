// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/io/pdb/DAminoAcidImportTests.cxxtest.hh
/// @brief  Tests that D-amino acids can be automatically identified by their chirality when the three-letter code is ambiguous.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>

// Protocols Headers (for conveneience, for setup):
#include <protocols/cyclic_peptide/DeclareBond.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("DAminoAcidImportTests");


class DAminoAcidImportTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
	}

	void tearDown(){

	}

	void test_import_ddpp_ldpp_cycpep(){
		TR << "Starting DAminoAcidImportTests:test_import_ddpp_ldpp_cycpep()." << std::endl;
		TR << "Test added 21 Dec 2020 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "Failure of this test means that Rosetta is not correctly assigning amino acid identities to residues that it imports (most likely to D-DPP or L-DPP, which have the same 3-letter codes, or to D-DAB, which has the same 3-letter code as L-DAB)." << std::endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/io/pdb/cycpep_with_ddpp_and_ldpp.pdb", core::import_pose::PDB_file );

		TS_ASSERT_EQUALS( pose.total_residue(), 10 );

		protocols::cyclic_peptide::DeclareBond decbond;
		decbond.set( 10, "C", 1, "N", false );
		decbond.apply(pose);

		TS_ASSERT_EQUALS( pose.residue_type(1).base_name(), "ASP" );
		TS_ASSERT_EQUALS( pose.residue_type(2).base_name(), "PRO" );
		TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "DASP" );
		TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "ASP" );
		TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "DDAB" );
		TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "DLEU" );
		TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "DASP" );
		TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "DPP" );
		TS_ASSERT_EQUALS( pose.residue_type(9).base_name(), "DDPP" );
		TS_ASSERT_EQUALS( pose.residue_type(10).base_name(), "DASP" );

		TR << "Completed DAminoAcidImportTests:test_import_ddpp_ldpp_cycpep()." << std::endl;
	}

	void test_import_ddpp_cycpep(){
		TR << "Starting DAminoAcidImportTests:test_import_ddpp_cycpep()." << std::endl;
		TR << "Test added 21 Dec 2020 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "Failure of this test means that Rosetta is not correctly assigning amino acid identities to residues that it imports (most likely to D-DPP, which has the same 3-letter code as L-DPP)." << std::endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/io/pdb/cycpep_with_ddpp.pdb", core::import_pose::PDB_file );

		TS_ASSERT_EQUALS( pose.total_residue(), 8 );

		protocols::cyclic_peptide::DeclareBond decbond;
		decbond.set( 8, "C", 1, "N", false );
		decbond.apply(pose);

		TS_ASSERT_EQUALS( pose.residue_type(1).base_name(), "PRO" );
		TS_ASSERT_EQUALS( pose.residue_type(2).base_name(), "AIB" );
		TS_ASSERT_EQUALS( pose.residue_type(3).base_name(), "GLU" );
		TS_ASSERT_EQUALS( pose.residue_type(4).base_name(), "DDPP" );
		TS_ASSERT_EQUALS( pose.residue_type(5).base_name(), "DPRO" );
		TS_ASSERT_EQUALS( pose.residue_type(6).base_name(), "LEU" );
		TS_ASSERT_EQUALS( pose.residue_type(7).base_name(), "DGLU" );
		TS_ASSERT_EQUALS( pose.residue_type(8).base_name(), "TRP" );

		TR << "Completed DAminoAcidImportTests:test_import_ddpp_cycpep()." << std::endl;
	}



};
