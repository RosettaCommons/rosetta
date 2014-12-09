// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/util/metalloprotein_import_test.cxxtest.hh
/// @brief  Test suite for import of metalloproteins using the -in:auto_setup_metals flag.
/// @details  The -in:auto_setup_metals flag is intended for automatic setup of covalent bonds and 
/// distance and angle constraints between metal-binding residues and metals in metalloprotein
/// structures.  Failure of this unit test means that something has been changed which prevents
/// this automatic setup.  Relevant files include:
///    -- src/core/util/metalloproteins_util.cc
///    -- src/core/util/metalloproteins_util.hh
///    -- src/core/conformation/Residue.cc
///    -- src/core/conformation/Residue.hh
///    -- src/core/conformation/ResidueType.cc
///    -- src/core/conformation/ResidueType.hh
/// @author Vikram K. Mulligan

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>


using namespace core;

class metalloprotein_import : public CxxTest::TestSuite
{

public:
	metalloprotein_import() {}

	// Shared initialization goes here.
	void setUp() {
		//core_init(); //Tests should FAIL without the -in:auto_setup_metals flag.
		core_init_with_additional_options( "-in:auto_setup_metals" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Tests that metalloproteins are imported properly, with
	/// automatic setup of covalent bonds.
	void test_metalloprotein_import() {
		pose::Pose pose;
		const std::string original_file_name("core/util/2c9v_stripped.pdb");
		import_pose::pose_from_pdb(pose, original_file_name);

		//Check that all the residues that should be bonded are bonded:
		TS_ASSERT( pose.residue(46).is_bonded(154) );
		TS_ASSERT( pose.residue(48).is_bonded(154) );
		TS_ASSERT( pose.residue(63).is_bonded(154) );
		TS_ASSERT( pose.residue(120).is_bonded(154) );

		TS_ASSERT( pose.residue(63).is_bonded(155) );
		TS_ASSERT( pose.residue(71).is_bonded(155) );
		TS_ASSERT( pose.residue(80).is_bonded(155) );
		TS_ASSERT( pose.residue(83).is_bonded(155) );

		TS_ASSERT( pose.residue(201).is_bonded(309) );
		TS_ASSERT( pose.residue(203).is_bonded(309) );
		TS_ASSERT( pose.residue(218).is_bonded(309) );
		TS_ASSERT( pose.residue(275).is_bonded(309) );

		TS_ASSERT( pose.residue(218).is_bonded(310) );
		TS_ASSERT( pose.residue(226).is_bonded(310) );
		TS_ASSERT( pose.residue(235).is_bonded(310) );
		TS_ASSERT( pose.residue(238).is_bonded(310) );

	}
};
