// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/AACompositionEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::AACompositionEnergy
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/methods/AACompositionEnergySetup.hh>
#include <core/scoring/methods/AACompositionEnergy.hh>

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
#include <core/conformation/util.hh>

#include <core/pose/annotated_sequence.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.scoring.methods.AACompositionEnergy.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;
using namespace core::scoring::annealing;

using namespace core::pack;
using namespace core::pack::task;
using namespace core::pack::rotamer_set;

class AACompositionEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation using the trp cage with noncanonicals.
	/// @details This test checks that we can impose the requirement that a pose contain exactly
	/// three trans-ACPC residues.
	void test_energy_eval_exactly_three_transACPC() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/methods/exactly_three_transACPC.comp -mute all -unmute core.scoring.methods.AACompositionEnergy.cxxtest");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests::test_energy_eval_exactly_three_transACPC()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose contain exactly three cisACPC." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Append one more transACPC:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("transACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd1, true, 0, 20, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+1ACPC:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Append one more transACPC:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("transACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd2, true, 0, 21, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 22 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 22 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+2ACPC:\t" << "20.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 20.0, 1e-6 );

		//Append one more transACPC:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 22 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 22 );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("transACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd3, true, 0, 22, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+3ACPC:\t" << "0.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.0, 1e-6 );

		//Append one more transACPC:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("transACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd4, true, 0, 23, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+4ACPC:\t" << "45.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 45.0, 1e-6 );

		//Append one more transACPC:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );
		core::conformation::ResidueOP new_rsd5( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("transACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd5, true, 0, 24, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 25 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 25 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+5ACPC:\t" << "45.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 45.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests::test_energy_eval_exactly_three_transACPC() complete." << std::endl;
			TR.flush();
		}
		return;
	}


	/// @brief Test the energy calculation using the trp cage.
	/// @details This test checks that we can impose the requirement that a pose contain exactly
	/// one tryptophan using this scoring term.
	void test_energy_eval_exactly_one_trp() {
		core_init_with_additional_options("-score:aa_composition_setup_file exactly_one_trp.comp -mute all -unmute core.scoring.methods.AACompositionEnergy.cxxtest");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests::test_energy_eval_exactly_one_trp()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose contain exactly one trp." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.0, 1e-6 );

		//Append one more tryptophan:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TRP") ) );
		trpcage.append_residue_by_bond(*new_rsd1, true, 0, 20, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+1trp:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Append one more tryptophan:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TRP") ) );
		trpcage.append_residue_by_bond(*new_rsd2, true, 0, 21, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 22 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 22 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+2trp:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Create another trp cage:
		Pose trpcage2( create_trpcage_ideal_pose() );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TYR") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 6 ), *new_rsd3, trpcage2.conformation(), true);
		trpcage2.replace_residue( 6, *new_rsd3, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-trp:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests::test_energy_eval_exactly_one_trp() complete." << std::endl;
			TR.flush();
		}
		return;
	}

	/// @brief Test the energy calculation using the trp cage.
	/// @details This test checks that we can impose the requirement that a pose contain exactly two aromatic residues.
	void test_energy_eval_two_aromatics() {
		core_init_with_additional_options("-score:aa_composition_setup_file two_aromatics.comp -mute all -unmute core.scoring.methods.AACompositionEnergy.cxxtest");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests::test_energy_eval_two_aromatics()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose have exactly two aromatics." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.0, 1e-6 );

		//Add one more tryptophan:
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TRP") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 7 ), *new_rsd1, trpcage.conformation(), true);
		trpcage.replace_residue( 7, *new_rsd1, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+trp:\t" << "60.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 60.0, 1e-6 );

		//Add one more phenylalanine:
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("PHE") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 8 ), *new_rsd2, trpcage.conformation(), true);
		trpcage.replace_residue( 8, *new_rsd2, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+trp+phe:\t" << "60.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 60.0, 1e-6 );

		//Create another trp cage and mutate the trp to ala:
		Pose trpcage2( create_trpcage_ideal_pose() );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 6 ), *new_rsd3, trpcage2.conformation(), true);
		trpcage2.replace_residue( 6, *new_rsd3, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-trp:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		//Mutate the tyr to ala:
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 3 ), *new_rsd4, trpcage2.conformation(), true);
		trpcage2.replace_residue( 3, *new_rsd4, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-trp-tyr:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests::test_energy_eval_two_aromatics() complete." << std::endl;
			TR.flush();
		}
		return;
	}

	/// @brief Test the energy calculation using the trp cage.
	/// @details This test checks that we can impose the requirement that a pose contain 10% aromatic residues.
	void test_energy_eval_ten_percent_aromatic() {
		core_init_with_additional_options("-score:aa_composition_setup_file ten_percent_aromatic.comp -mute all -unmute core.scoring.methods.AACompositionEnergy.cxxtest");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests::test_energy_eval_ten_percent_aromatic()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose be 10 percent aromatics." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.0, 1e-6 );

		//Add one more tryptophan:
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TRP") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 7 ), *new_rsd1, trpcage.conformation(), true);
		trpcage.replace_residue( 7, *new_rsd1, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+trp:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Add one more phenylalanine:
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("PHE") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 8 ), *new_rsd2, trpcage.conformation(), true);
		trpcage.replace_residue( 8, *new_rsd2, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+trp+phe:\t" << "50.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.0, 1e-6 );

		//Create another trp cage and mutate the trp to ala:
		Pose trpcage2( create_trpcage_ideal_pose() );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 6 ), *new_rsd3, trpcage2.conformation(), true);
		trpcage2.replace_residue( 6, *new_rsd3, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-trp:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		//Mutate the tyr to ala:
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 3 ), *new_rsd4, trpcage2.conformation(), true);
		trpcage2.replace_residue( 3, *new_rsd4, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-trp-tyr:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests::test_energy_eval_ten_percent_aromatic() complete." << std::endl;
			TR.flush();
		}
		return;
	}

	/// @brief Test the energy calculation using the trp cage.
	/// @details This test checks that we can impose the requirement that a pose contain 20% proline.
	void test_energy_eval_twenty_percent_pro() {
		core_init_with_additional_options("-score:aa_composition_setup_file twenty_percent_pro.comp -mute all -unmute core.scoring.methods.AACompositionEnergy.cxxtest");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests::test_energy_eval_twenty_percent_pro()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose be 20 percent prolines." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.0, 1e-6 );

		//Add one more proline:
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("PRO") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 2 ), *new_rsd1, trpcage.conformation(), true);
		trpcage.replace_residue( 2, *new_rsd1, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+pro:\t" << "10.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 10.0, 1e-6 );

		//Add one more proline:
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("PRO") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 3 ), *new_rsd2, trpcage.conformation(), true);
		trpcage.replace_residue( 3, *new_rsd2, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+2pro:\t" << "20.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 20.0, 1e-6 );

		//Add one more proline:
		core::conformation::ResidueOP new_rsd2b( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("PRO") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage.residue( 4 ), *new_rsd2b, trpcage.conformation(), true);
		trpcage.replace_residue( 4, *new_rsd2b, false );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+3pro:\t" << "20.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 20.0, 1e-6 );

		//Create another trp cage and mutate a pro to ala:
		Pose trpcage2( create_trpcage_ideal_pose() );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 17 ), *new_rsd3, trpcage2.conformation(), true);
		trpcage2.replace_residue( 17, *new_rsd3, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-pro:\t" << "25.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 25.0, 1e-6 );

		//Mutate the another pro to ala:
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 18 ), *new_rsd4, trpcage2.conformation(), true);
		trpcage2.replace_residue( 18, *new_rsd4, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-2pro:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		//Mutate the another pro to ala:
		core::conformation::ResidueOP new_rsd5( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 19 ), *new_rsd5, trpcage2.conformation(), true);
		trpcage2.replace_residue( 19, *new_rsd5, false );

		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-3pro:\t" << "50.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 50.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests::test_energy_eval_twenty_percent_pro() complete." << std::endl;
			TR.flush();
		}
		return;
	}

	/// @brief Test the energy calculation with the packer.
	/// @author Alex Ford.
	void test_energy_annealing( ) {
		core_init_with_additional_options("-score:aa_composition_setup_file exactly_one_ala.comp");

		// Setup score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		// Setup test pose
		Pose pose;
		make_pose_from_sequence( pose, "AGGGGGGG", "fa_standard");
		TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);

		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");

		PackerEnergy prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 100, 1e-6);

		// Setup packer task and packer objects
		PackerTaskOP task( TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;

		for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		TS_ASSERT_EQUALS( rotsets->nrotamers(), pose.total_residue() * 2);

		core::pack::interaction_graph::ResidueArrayAnnealingEvaluator ev;
		ev.initialize( scorefxn, pose, *rotsets, packer_neighbor_graph);

		TS_ASSERT_EQUALS( ev.get_num_nodes(), 8);
		TS_ASSERT_EQUALS( ev.get_num_total_states(), 16);
		TS_ASSERT_EQUALS( ev.get_num_states_for_node(1), 2);

		// Base assignment should be base pose score...
		TS_ASSERT_EQUALS( ev.any_vertex_state_unassigned(), true );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		// Test state assignment and consideration
		ev.set_state_for_node(1, 1);
		for ( int r = 2; r <= ev.get_num_nodes(); ++r ) {
			ev.set_state_for_node( r, 2);
		}

		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.set_state_for_node(2, 1);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		PackerEnergy delta_energy;
		PackerEnergy pre_energy;

		ev.consider_substitution(2, 2, delta_energy, pre_energy);

		TS_ASSERT_DELTA( delta_energy, -100, 1e-6);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		ev.commit_considered_substitution();
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		// Test via pack_rotamers run
		pack_rotamers(pose, scorefxn, task);
		PackerEnergy postpack_score = scorefxn(pose);
		TS_ASSERT_DELTA(postpack_score, 0, 1e-6);

		std::map< std::string, int > aa_count;
		aa_count["ALA"] = 0;
		aa_count["GLY"] = 0;

		for ( core::Size r = 1; r<= pose.total_residue(); ++r ) {
			aa_count[ pose.residue(r).name3() ] += 1;
		}

		TS_ASSERT_EQUALS( aa_count["ALA"], 1);
		TS_ASSERT_EQUALS( aa_count["GLY"], pose.total_residue() - 1);
	}

	/// @brief Test the tail functions with constant below, linear above.
	///
	void test_tailfunctions_const_lin() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/methods/tailfunction_const.comp");

		if ( TR.visible() ) {
			TR << "Starting test_tailfunctions_const_lin()." << std::endl;
		}

		// Set up score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		{ //Scope 1: within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);
		}

		{ //Scope 2: one below ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 3: two below ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 4: three below ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAGGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 5: one above ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: two above ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 40, 1e-6);
		}

		{ //Scope 6: three above ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAAGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 60, 1e-6);
		}

		if ( TR.visible() ) {
			TR << "Test test_tailfunctions_const_lin() complete." << std::endl;
			TR.flush();
		}

		return;
	}

	/// @brief Test the tail functions with linear below, const above.
	///
	void test_tailfunctions_lin_const() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/methods/tailfunction_linear.comp");

		if ( TR.visible() ) {
			TR << "Starting test_tailfunctions_lin_const()." << std::endl;
		}

		// Set up score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		{ //Scope 1: within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);
		}

		{ //Scope 2: one below ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 3: two below ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 4: three below ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAGGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 30, 1e-6);
		}

		{ //Scope 5: one above ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: two above ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: three above ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAAGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		if ( TR.visible() ) {
			TR << "Test test_tailfunctions_lin_const() complete." << std::endl;
			TR.flush();
		}

		return;
	}

	/// @brief Test the tail functions with quadratic above and below.
	///
	void test_tailfunctions_quadratic() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/methods/tailfunction_quadratic.comp");

		if ( TR.visible() ) {
			TR << "Starting test_tailfunctions_quadratic()." << std::endl;
		}

		// Set up score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		{ //Scope 1: within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);
		}

		{ //Scope 2: one below ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 3: two below ideal, out of test range (quadratic region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 40, 1e-6);
		}

		{ //Scope 4: three below ideal, out of test range (quadratic region)
			Pose pose;
			make_pose_from_sequence( pose, "AAGGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 90, 1e-6);
		}

		{ //Scope 5: one above ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: two above ideal, out of test range (quadratic region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 80, 1e-6);
		}

		{ //Scope 6: three above ideal, out of test range (quadratic region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAAGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 180, 1e-6);
		}

		if ( TR.visible() ) {
			TR << "Test test_tailfunctions_quadratic() complete." << std::endl;
			TR.flush();
		}

		return;
	}

};
