
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/aa_composition/MHCEpitopeEnergy.cxxtest.hh
/// @brief  Test suite for protocols::aa_composition::::AddMHCEptiopeConstraintMover, which adds MHC constraints to
/// @brief  specific parts of a pose.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <protocols/aa_composition/AddMHCEpitopeConstraintMover.hh>
#include <protocols/aa_composition/ClearCompositionConstraintsMover.hh>

#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

// Package Headers
#include <test/core/init_util.hh>

#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("protocols.aa_composition.AddMHCEpitopeConstraintMover.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace protocols::aa_composition;
using namespace core::scoring;
using namespace core::pose;

class AddMHCEpitopeConstraintMoverTests : public CxxTest::TestSuite {

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

	/// @brief Test the AddMHCEpitopeConstraintMover and the ClearCompositionConstraintsMover
	/// @author Brahm Yachnin
	void test_mhc_csts() {
		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");

		//Associate the config with the scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );

		//Calculate the energy of the starting pose
		//Should be 0, since there are no csts and no generalized configuration
		TS_ASSERT_EQUALS( scorefxn(pose), 0 );

		//Create an AddMHCEpitopeConstraintMover object
		AddMHCEpitopeConstraintMoverOP add_mhc_csts( utility::pointer::make_shared<AddMHCEpitopeConstraintMover>() );
		//Create the csts from file contents
		std::string config_string = "method matrix propred8\nalleles *\nthreshold 5";
		add_mhc_csts->create_constraint_from_file_contents( config_string );
		//Create a residue selector for the first 9 residues (one allele).
		core::select::residue_selector::ResidueIndexSelectorCOP selector( utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>("1-9") );

		//Associate that residue selector with the mover.
		add_mhc_csts->add_residue_selector( selector );

		//Apply add_mhc_csts to the pose
		add_mhc_csts->apply(pose);

		core::Real base_cst_score = scorefxn(pose);
		TS_ASSERT_EQUALS( base_cst_score, 8.0 );

		//Create the clear csts mover
		ClearCompositionConstraintsMoverOP clear_mhc_csts( utility::pointer::make_shared<ClearCompositionConstraintsMover>() );
		clear_mhc_csts->apply(pose);
		TS_ASSERT_EQUALS( scorefxn(pose), 0 );

		//Create another AddMHCEpitopeConstraintMover object, as before
		AddMHCEpitopeConstraintMoverOP add_weighted_mhc_csts( utility::pointer::make_shared<AddMHCEpitopeConstraintMover>() );
		add_weighted_mhc_csts->create_constraint_from_file_contents( config_string );
		add_weighted_mhc_csts->add_residue_selector( selector );
		//For this cst, apply a weight of 0.5
		core::Real weight = 0.5;
		add_weighted_mhc_csts->add_weight(weight);
		add_weighted_mhc_csts->apply(pose);

		//The score should now be weight * the base score
		TS_ASSERT_EQUALS( scorefxn(pose), base_cst_score * weight );

		//Clear the csts again
		clear_mhc_csts->apply(pose);
		TS_ASSERT_EQUALS( scorefxn(pose), 0);

		//Now re-apply both csts at the same time
		add_mhc_csts->apply(pose);
		add_weighted_mhc_csts->apply(pose);

		//The score should now be (1+weight) * the base score
		TS_ASSERT_EQUALS( scorefxn(pose), base_cst_score * (weight+1) );

		//Clear the csts again
		clear_mhc_csts->apply(pose);
		TS_ASSERT_EQUALS( scorefxn(pose), 0);

		TR << "End of test_mhc_csts." << std::endl;
	}

};
