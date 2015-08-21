// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/AARepeatEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::AARepeatEnergy
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/methods/AARepeatEnergy.hh>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.scoring.methods.AARepeatEnergy.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class AARepeatEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation for varying numbers of repeating residues.
	///
	void test_energy_eval() {
		if ( TR.visible() ) {
			TR << "Starting AARepeatEnergyTests::test_energy_eval()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_repeat_energy score term evaluates its energy correctly." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_repeat_energy, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.5, 1e-6 );

		//Append three more alanines:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd1, true, 0, 20, 0, false);
		trpcage.append_residue_by_bond(*new_rsd2, true, 0, 21, 0, false);
		trpcage.append_residue_by_bond(*new_rsd3, true, 0, 22, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+3:\t" << "1.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 1.0, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd4, true, 0, 23, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+4:\t" << "5.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 5.5, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );
		core::conformation::ResidueOP new_rsd5( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd5, true, 0, 24, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 25 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 25 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+5:\t" << "50.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.5, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 25 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 25 );
		core::conformation::ResidueOP new_rsd6( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd6, true, 0, 25, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 26 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 26 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+6:\t" << "50.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.5, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AARepeatEnergyTests::test_energy_eval() complete." << std::endl;
			TR.flush();
		}
		return;
	}

};


