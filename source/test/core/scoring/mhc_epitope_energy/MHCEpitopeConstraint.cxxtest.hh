// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEpitopeConstraint
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pose/Pose.hh>

// Package Headers
#include <test/core/init_util.hh>

#include <core/pose/annotated_sequence.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopeConstraint.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core::pose;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopeConstraintTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //
	// Shared variables
	Pose pose;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		// Setup test pose
		std::string sequence = "YFCTRAFRILAWIGIQNPTS";
		make_pose_from_sequence( pose, sequence, "fa_standard");
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Some unit tests to ensure that the MHCEpitopeConstraint class can be setup correctly.
	/// @author Brahm Yachnin
	void test_mhc_constraints() {
		// Set some parameters to use
		core::Real cstweight = 0.5;

		// Create a residue selector for the first 9 residues.
		core::select::residue_selector::ResidueIndexSelectorCOP selector( utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>("1-9") );

		// MHCEpitopeSetup string
		std::string setup_string = "method matrix propred8\nalleles *\nthreshold 5";
		// Make our own MHCEpitopeEnergySetup object (not a cst)
		MHCEpitopeEnergySetupOP mhc_setup_obj( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		mhc_setup_obj->initialize_from_file_contents(setup_string);

		// Make a new MHCEpitopeConstraint object
		MHCEpitopeConstraintOP mhc_cst_obj( utility::pointer::make_shared<MHCEpitopeConstraint>() );

		// Set the desired parameters
		mhc_cst_obj->initialize_from_file_contents( setup_string );
		mhc_cst_obj->set_selector( selector );
		mhc_cst_obj->set_cst_weight( cstweight );

		// Check that the parameters all match up appropriately
		TS_ASSERT_EQUALS( mhc_cst_obj->selector()->apply(pose), selector->apply(pose) );
		TS_ASSERT_EQUALS( mhc_cst_obj->mhc_epitope_energy_setup()->report(), mhc_setup_obj->report() );
		TS_ASSERT_EQUALS( mhc_cst_obj->get_cst_weight(), cstweight);

		TR << "End of test_mhc_constraints." << std::endl;
	}
};
