// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/component_reading.cxxtest.hh
/// @brief  test suite for the functionality of reading ResidueTypes from components file.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Program headers

// Basic headers
#include <basic/Tracer.hh>

// Utility headers

// ObjexxFCL headers

// C++ headers
#include <string>

static basic::Tracer TR("core.io.pdb.component_reading.cxxtest");

class component_reading_Tests : public CxxTest::TestSuite
{

public:
	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-ignore_waters -PDB_components_file core/chemical/mmCIF/components_trimmed.cif" );
	}

	/// @brief Test the scoring scheme
	void test_load_components() {
		core::pose::PoseOP pose;

		pose = core::import_pose::pose_from_pdb("core/import_pose/3KBA.pdb");

		TS_ASSERT_EQUALS( pose->total_residue(), 496 ); // No Waters
		TS_ASSERT_EQUALS( pose->residue(244).name3(), "WOW"); // Actually the 494th residue in file
		TS_ASSERT_EQUALS( pose->residue(245).name3(), "SO4"); // Acutally the 495th residue in file
		TS_ASSERT_EQUALS( pose->residue(496).name3(), "WOW");
		TS_ASSERT_EQUALS( pose->residue(244).natoms(), 49 ); // WOW
		TS_ASSERT_EQUALS( pose->residue(244).nheavyatoms(), 27 ); // WOW


	}

};

