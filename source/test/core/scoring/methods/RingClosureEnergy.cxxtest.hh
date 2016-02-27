// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/RingClosureEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::RingClosureEnergy.cc
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/RingClosureEnergy.hh>
#include <platform/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>


// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>

#include <core/kinematics/DomainMap.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//C++ Headers
#include <math.h>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::id;
using namespace core::kinematics;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

static basic::Tracer TR("core.scoring.methods.RingClosureEnergy.cxxtest");

class RingClosureEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using namespace core;
		core_init_with_additional_options( "-write_all_connect_info -connect_info_cutoff 0 -output_virtual true" );
	}

	// Shared finalization goes here.
	void tearDown() {}

	// --------------- Test Cases --------------- //
	void test_eval_energy_proline()
	{
		if ( TR.visible() ) {
			TR << "Starting RingClosureEnergyTests::test_eval_energy_proline()." << std::endl;
			TR << "This test applies the ring_close energy term to the trp cage, then confirms that the energy is zero for non-proline residues and equal to expected values for proline values.  Expected values are based on the intra-residue part of the pro_close energy." << std::endl;
			TR << "Test created by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory, 11 July 2015." << std::endl;
		}

		Pose trpcage( create_trpcage_ideal_pose() );

		RingClosureEnergy ringclose_energy;

		utility::vector1 <Real> expected_energies;
		for ( core::Size i=1; i<=11; ++i ) expected_energies.push_back(0.0);
		expected_energies.push_back(0.0124967778684113); //12
		for ( core::Size i=13; i<=16; ++i ) expected_energies.push_back(0.0);
		expected_energies.push_back(0.2063039233793921); //17
		expected_energies.push_back(0.1033872216433958); //18
		expected_energies.push_back(0.1980550619309646); //19
		expected_energies.push_back(0.0); //20

		EnergyMap emap;
		for ( core::Size i=1; i<=20; ++i ) {
			emap.zero();
			ringclose_energy.residue_energy( trpcage.residue( i ), trpcage, emap );
			if ( TR.visible() ) TR << "Ring_close energy of residue " << i << ": " << emap[ring_close] << std::endl;
			TS_ASSERT_DELTA( emap[ ring_close ], expected_energies[i], 1e-12 );
		}

		if ( TR.visible() ) {
			TR << "Finished RingClosureEnergyTests::test_eval_energy_proline()." << std::endl;
			TR.flush();
		}

		return;
	}

	void test_ringclose_deriv_check_w_total_flexibility()
	{
		if ( TR.visible() ) {
			TR << "Starting RingClosureEnergyTests::test_ringclose_deriv_check_w_total_flexibility()." << std::endl;
			TR << "This test checks that the RingClosureEnergy derivatives, computed analytically by the function's internal derivative calculation functions, match the numerical derivatives." << std::endl;
			TR << "Test created by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory, 11 July 2015." << std::endl;
		}
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( ring_close, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.simple_deriv_check( false, 1e-6 );

		if ( TR.visible() ) {
			TR << "Finished RingClosureEnergyTests::test_ringclose_deriv_check_w_total_flexibility()." << std::endl;
			TR.flush();
		}
	}

	void test_eval_energy_cisACPC()
	{
		if ( TR.visible() ) {
			TR << "Starting RingClosureEnergyTests::test_eval_energy_cisACPC()." << std::endl;
			TR << "This test adds a cisACPC residue to the trp cage, applies the ring_close energy term, then confirms that the energy equal to expected values for cisACPC.  Expected values are based on manual calculation." << std::endl;
			TR << "Test created by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory, 11 July 2015." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );


		Pose trpcage( create_trpcage_ideal_pose() );

		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );

		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("cisACPC") ) );
		trpcage.append_residue_by_bond(*new_rsd, true, 0, 20, 0, false);

		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );

		EnergyMap emap;
		RingClosureEnergy ringclose_energy;

		{ //Scope1: initial conformation
			//trpcage.dump_pdb("vtemp1.pdb"); //DELETE ME

			Distance const d1( trpcage.residue(21).xyz("VCM").distance( trpcage.residue(21).xyz("CM") ) );
			Distance const d2( trpcage.residue(21).xyz("VCD").distance( trpcage.residue(21).xyz("CD") ) );
			Real const predicted_energy(100*(pow(d1,2.0)+pow(d2,2.0)));

			if ( TR.visible() ) {
				TR << "VCM-CM distance:\t" << d1 << "\t" << pow(d1,2.0) << std::endl;
				TR << "VCD-CD distance:\t" << d2 << "\t" << pow(d2,2.0) << std::endl;
				TR << "Predicted ring closure energy:\t" << predicted_energy << std::endl;
			}

			emap.zero();
			ringclose_energy.residue_energy( trpcage.residue( 21 ), trpcage, emap );
			if ( TR.visible() ) TR << "Ring_close energy of residue 21 (cisACPC): " << emap[ring_close] << std::endl;
			TS_ASSERT_DELTA( predicted_energy, emap[ring_close], 1e-6 );
		}

		{ //Scope2: altered conformation
			if ( TR.visible() ) TR << "PERTURBING TRP CAGE POSE" << std::endl;
			trpcage.set_chi(1,21,53.7);  //Arbitrarily chosen value
			trpcage.set_chi(2,21,-143.2);  //Arbitrarily chosen value
			//trpcage.dump_pdb("vtemp2.pdb"); //DELETE ME

			Distance const d1( trpcage.residue(21).xyz("VCM").distance( trpcage.residue(21).xyz("CM") ) );
			Distance const d2( trpcage.residue(21).xyz("VCD").distance( trpcage.residue(21).xyz("CD") ) );
			Real const predicted_energy(100*(pow(d1,2.0)+pow(d2,2.0)));

			if ( TR.visible() ) {
				TR << "VCM-CM distance:\t" << d1 << "\t" << pow(d1,2.0) << std::endl;
				TR << "VCD-CD distance:\t" << d2 << "\t" << pow(d2,2.0) << std::endl;
				TR << "Predicted ring closure energy:\t" << predicted_energy << std::endl;
			}

			emap.zero();
			ringclose_energy.residue_energy( trpcage.residue( 21 ), trpcage, emap );
			if ( TR.visible() ) TR << "Ring_close energy of residue 21 (cisACPC): " << emap[ring_close] << std::endl;
			TS_ASSERT_DELTA( predicted_energy, emap[ring_close], 1e-6 );
		}

		if ( TR.visible() ) {
			TR << "Finished RingClosureEnergyTests::test_eval_energy_cisACPC()." << std::endl;
			TR.flush();
		}

		return;
	}

	void test_eval_energy_cisACHC()
	{
		if ( TR.visible() ) {
			TR << "Starting RingClosureEnergyTests::test_eval_energy_cisACHC()." << std::endl;
			TR << "This test adds a cisACHC residue to the trp cage, applies the ring_close energy term, then confirms that the energy equal to expected values for cisACHC.  Expected values are based on manual calculation." << std::endl;
			TR << "Test created by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory, 11 July 2015." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );


		Pose trpcage( create_trpcage_ideal_pose() );

		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );

		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("cisACHC") ) );
		trpcage.append_residue_by_bond(*new_rsd, true, 0, 20, 0, false);

		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 21 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 21 );

		EnergyMap emap;
		RingClosureEnergy ringclose_energy;

		{ //Scope1: initial conformation
			//trpcage.dump_pdb("vtemp1.pdb"); //DELETE ME

			Distance const d1( trpcage.residue(21).xyz("VCM").distance( trpcage.residue(21).xyz("CM") ) );
			Distance const d2( trpcage.residue(21).xyz("VCE").distance( trpcage.residue(21).xyz("CE") ) );
			Real const predicted_energy(100*(pow(d1,2.0)+pow(d2,2.0)));

			if ( TR.visible() ) {
				TR << "VCM-CM distance:\t" << d1 << "\t" << pow(d1,2.0) << std::endl;
				TR << "VCE-CE distance:\t" << d2 << "\t" << pow(d2,2.0) << std::endl;
				TR << "Predicted ring closure energy:\t" << predicted_energy << std::endl;
			}

			emap.zero();
			ringclose_energy.residue_energy( trpcage.residue( 21 ), trpcage, emap );
			if ( TR.visible() ) TR << "Ring_close energy of residue 21 (cisACHC): " << emap[ring_close] << std::endl;
			TS_ASSERT_DELTA( predicted_energy, emap[ring_close], 1e-6 );
		}

		{ //Scope2: altered conformation
			if ( TR.visible() ) TR << "PERTURBING TRP CAGE POSE" << std::endl;
			trpcage.set_chi(1,21,53.7);  //Arbitrarily chosen value
			trpcage.set_chi(2,21,-143.2);  //Arbitrarily chosen value
			//trpcage.dump_pdb("vtemp2.pdb"); //DELETE ME

			Distance const d1( trpcage.residue(21).xyz("VCM").distance( trpcage.residue(21).xyz("CM") ) );
			Distance const d2( trpcage.residue(21).xyz("VCE").distance( trpcage.residue(21).xyz("CE") ) );
			Real const predicted_energy(100*(pow(d1,2.0)+pow(d2,2.0)));

			if ( TR.visible() ) {
				TR << "VCM-CM distance:\t" << d1 << "\t" << pow(d1,2.0) << std::endl;
				TR << "VCE-CE distance:\t" << d2 << "\t" << pow(d2,2.0) << std::endl;
				TR << "Predicted ring closure energy:\t" << predicted_energy << std::endl;
			}

			emap.zero();
			ringclose_energy.residue_energy( trpcage.residue( 21 ), trpcage, emap );
			if ( TR.visible() ) TR << "Ring_close energy of residue 21 (cisACHC): " << emap[ring_close] << std::endl;
			TS_ASSERT_DELTA( predicted_energy, emap[ring_close], 1e-6 );
		}

		if ( TR.visible() ) {
			TR << "Finished RingClosureEnergyTests::test_eval_energy_cisACHC()." << std::endl;
			TR.flush();
		}

		return;
	}


};


