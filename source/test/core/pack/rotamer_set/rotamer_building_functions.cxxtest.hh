// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 rotamer_building_functions.cxxtest.hh
/// @brief   Test suite for rotamer building functions
/// @author  Labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/pack/rotamer_set/rotamer_building_functions.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>


// Utility header
#include <utility/vector1.hh>

// C++ header
#include <string>


class RotamerBuildingFunctionsTests : public CxxTest::TestSuite {
public:
	// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


	// Tests //////////////////////////////////////////////////////////////////////////////////////////////////////////
	void test_build_rotamers_from_rotamer_bins()
	{
		using namespace utility;
		using namespace core::conformation;
		using namespace core::pack::rotamer_set;
		using namespace std;
		using namespace utility;
		using namespace core::chemical;
		using namespace core::chemical::orbitals;
		using namespace core::pose;
		using namespace core::import_pose;

		// Load params file for glycerol to use as a generic test case.
		vector1<string> params;
		params.push_back("core/pack/rotamer_set/GOL.params");

 		// Instantiate default chemical type sets.
		ResidueTypeSet & res_type_set = ChemicalManager::get_instance()->nonconst_residue_type_set("fa_standard");

		// Generate a non-standard ResidueTypeSet that includes parameters for glycerol.
		res_type_set.read_files(params);

		// Load a pose with a single glycerol residue.
		Pose glycerol;
		pose_from_pdb(glycerol, res_type_set, "core/pack/rotamer_set/glycerol.pdb");

		// Extract out the residue for rotamer testing.
		Residue const & GOL = glycerol.residue(1);  // C3H8O3

		vector1<ResidueOP> rotamers;

		build_rotamers_from_rotamer_bins(GOL, rotamers);

		// Glycerol has 3 hydroxyl groups and, thus, 3 chi angles.  The params file specifies 3 rotamer bins for each
		// chi.  Thus, the number of romaters generated should be 3^3 = 27.
		TS_ASSERT_EQUALS(rotamers.size(), 27);

		// TODO: Test for further cases, such as when a chi angle has no rotamer bins specified.
	}
};  // class RotamerBuildingFunctionsTests
