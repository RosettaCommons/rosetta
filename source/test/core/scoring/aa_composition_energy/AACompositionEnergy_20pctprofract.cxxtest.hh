// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/aa_composition_energy/AACompositionEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::aa_composition_energy::AACompositionEnergy, an energy term for controlling
/// sequence composition during design.
/// @details See also the core::conformation::symmetry::MirrorSymmetricConformation unit tests.  These have
/// another example of AAComposition being set up from code (with constraints attached to the pose).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergy.hh>

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


static basic::Tracer TR("core.scoring.aa_composition_energy.AACompositionEnergy.cxxtest");

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

class AACompositionEnergyTests_20pctprofract : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation using the trp cage.
	/// @details This test checks that we can impose the requirement that a pose contain 20% proline, with penalty ranges specified as fractions.
	void test_energy_eval_twenty_percent_pro_fract_ranges() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/aa_composition_energy/twenty_percent_pro_fract_ranges.comp -out:levels core.scoring.aa_composition_energy.AACompositionEnergy:500");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests_20pctprofract::test_energy_eval_twenty_percent_pro_fract_ranges()." << std::endl;
			TR << "Test created 28 Apr. 2016 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that requires that a pose be 20 percent prolines.  This test uses a composition file that specifies ranges as fractions." << std::endl;
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
		if ( TR.visible() ) TR << "TrpCage+3pro:\t" << "36.67\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 36.66666666666667, 1e-6 );

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
		if ( TR.visible() ) TR << "TrpCage-3pro:\t" << "75.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 75.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests_20pctprofract::test_energy_eval_twenty_percent_pro() complete." << std::endl;
			TR.flush();
		}
		return;
	}

};
