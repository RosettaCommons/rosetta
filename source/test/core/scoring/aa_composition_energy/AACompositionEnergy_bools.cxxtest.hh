// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

class AACompositionEnergyTests_bools : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation using the trp cage with a .comp file with fairly complex Boolean logic.
	/// @details This test defines a count group in which a residue is counted if it is a tryptophan OR it is ((charged or aliphatic) and not (negatively charged or argenine or leucine)).
	/// So the following residue types should be counted: AIKMPVW.
	void test_energy_eval_complex_boolean_logic() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/aa_composition_energy/complex_booleans.comp -out:levels core.scoring.aa_composition_energy.AACompositionEnergy:500");
		if ( TR.visible() ) {
			TR << "Starting AACompositionEnergyTests_bools::test_energy_eval_complex_boolean_logic()." << std::endl;
			TR << "Test created 21 Nov 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_composition score term evaluates its energy correctly.  It uses the trp cage, and scores using a setup file that counts a residue if it is a tryptophan OR it is ((charged or aliphatic) and not (negatively charged or argenine or leucine))." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		TR << "Trp cage sequence: " << trpcage.sequence() << std::endl;;
		ScoreFunction sfxn;
		sfxn.set_weight( aa_composition, 1.0 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "7.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 7.0, 1e-6 );

		//Mutate the trp to leu:
		Pose trpcage2(trpcage);
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("LEU") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 6 ), *new_rsd1, trpcage2.conformation(), true);
		trpcage2.replace_residue( 6, *new_rsd1, false );
		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-W6L:\t" << "6.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 6.0, 1e-6 );
		//Mutate a pro to glu:
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("GLU") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 19 ), *new_rsd2, trpcage2.conformation(), true);
		trpcage2.replace_residue( 19, *new_rsd2, false );
		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-W6L,P19E:\t" << "5.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 5.0, 1e-6 );
		//Add an arginine:
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ARG") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 4 ), *new_rsd3, trpcage2.conformation(), true);
		trpcage2.replace_residue( 4, *new_rsd3, false );
		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-W6L,P19E,I4R:\t" << "4.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 4.0, 1e-6 );
		//Add back a tryptophan:
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("TRP") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 1 ), *new_rsd4, trpcage2.conformation(), true);
		trpcage2.replace_residue( 1, *new_rsd4, false );
		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-N1W,W6L,P19E,I4R:\t" << "5.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 5.0, 1e-6 );
		//Ged rid of 2 prolines:
		core::conformation::ResidueOP new_rsd5( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("LEU") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 17 ), *new_rsd5, trpcage2.conformation(), true);
		trpcage2.replace_residue( 17, *new_rsd5, false );
		core::conformation::ResidueOP new_rsd6( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ASN") ) );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( trpcage2.residue( 18 ), *new_rsd6, trpcage2.conformation(), true);
		trpcage2.replace_residue( 18, *new_rsd6, false );
		sfxn(trpcage2);
		if ( TR.visible() ) TR << "TrpCage-N1W,W6L,P19E,I4R,P17L,P18N:\t" << "3.0\t" << trpcage2.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage2.energies().total_energy(), 3.0, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AACompositionEnergyTests_bools::test_energy_eval_complex_boolean_logic() complete." << std::endl;
			TR.flush();
		}
		return;
	}

};
