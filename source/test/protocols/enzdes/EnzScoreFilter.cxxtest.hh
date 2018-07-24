// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/enzdes/EnzScoreFilter.cxxtest.hh
/// @brief
/// @author Brahm Yachnin (brahm.yachnin@rutgers.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tracer header
#include <basic/Tracer.hh>

// Unit headers
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/enzdes/EnzFilters.hh>
//#include <protocols/enzdes/AddorRemoveCsts.cc>

#include <core/scoring/ScoreType.hh>
//#include <core/scoring/ScoreFunction.hh>

//#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <utility/vector1.hh>

// --------------- Test Class --------------- //

// using declarations
using namespace core;

static basic::Tracer TR( "protocols.enzdes.EnzScoreFilter" );

class EnzScoreFilterTests : public CxxTest::TestSuite {

public:

	Real TOLERATED_ERROR;
	pose::Pose pose;
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP cst_io;

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		using namespace core::chemical;
		using namespace protocols::enzdes;

		protocols_init_with_additional_options("-run:preserve_header -extra_res_fa protocols/enzdes/D2N.params");
		TOLERATED_ERROR = 0.1;

		//Read in the pose
		core::import_pose::pose_from_file( pose, "protocols/enzdes/ligtest_it.pdb", core::import_pose::PDB_file);

		//Read in the enzdes constraint file into cst_io.
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		cst_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO(const_residue_set) );
		cst_io->read_enzyme_cstfile("protocols/enzdes/ligtest_it.cst");
	}

	void tearDown() {}

	// --------------- Test Cases --------------- //

	// Unit test for calculation of whole_pose cst energy.
	void test_whole_pose_enzdes_csts() {
		using namespace protocols::enzdes;
		using namespace core::scoring;

		//Set up a scorefunction with csts
		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->reset();
		scorefxn->set_weight( scoring::atom_pair_constraint, 1.0);
		scorefxn->set_weight( scoring::angle_constraint, 1.0);
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0);

		//Set up variables containing the defaults for each EnzScoreFilter option.
		std::string default_resnum = EnzScoreFilter::default_value_for_resnum();
		std::string default_cstid = EnzScoreFilter::default_value_for_cstid();
		core::Real default_threshold = EnzScoreFilter::default_value_for_threshold();

		//Add the constraints to the pose
		cst_io->add_constraints_to_pose( pose, scorefxn, true );

		//Create an EnzScoreFilter object.
		//The first "true" argument sets whole_pose to true.
		//The second "true" argument sets is_cstE to true (which overrides the score_type setting).
		//fa_rep is a arbitrary scoretype that will be overridden by is_cstE=true.
		EnzScoreFilter wholepose_test( default_resnum, default_cstid, scorefxn, fa_rep, default_threshold, true, true );

		//Compute the enzdes-style energy for the whole pose
		core::Real compute_result = wholepose_test.compute(pose);
		core::Real expect_result = 80.0313;
		TR << "test_whole_pose_enzdes_csts has a result of " << compute_result << ", compared to expected result of " << expect_result << std::endl;
		TS_ASSERT_DELTA( compute_result, expect_result, TOLERATED_ERROR );
		TR << "End of test_whole_pose_enzdes_csts." << std::endl;
	}

	// Unit test for calculation of cst energy using a resnum specifier or cstid.
	void test_resnum_enzdes_csts() {
		using namespace protocols::enzdes;
		using namespace core::scoring;

		//Set up a scorefunction with csts
		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->reset();
		scorefxn->set_weight( scoring::atom_pair_constraint, 1.0);
		scorefxn->set_weight( scoring::angle_constraint, 1.0);
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0);

		//Set up variables containing the defaults for each EnzScoreFilter option.
		std::string default_resnum = EnzScoreFilter::default_value_for_resnum();
		std::string default_cstid = EnzScoreFilter::default_value_for_cstid();
		core::scoring::ScoreType default_scoretype = core::scoring::score_type_from_name( EnzScoreFilter::default_value_for_score_type() );
		core::Real default_threshold = EnzScoreFilter::default_value_for_threshold();
		bool default_wholepose = EnzScoreFilter::default_value_for_whole_pose();

		//Add the constraints to the pose
		cst_io->add_constraints_to_pose( pose, scorefxn, true );

		//Create an EnzScoreFilter object.
		//The first argument specifies the desired resnum to use.
		//The "true" argument sets is_cstE to true (which overrides the score_type setting).
		EnzScoreFilter pdbnum_test( "1X", default_cstid, scorefxn, default_scoretype, default_threshold, default_wholepose, true );

		//Compute the enzdes-style energy for the given pdbnum
		core::Real compute_result = pdbnum_test.compute(pose);
		core::Real expect_result = 40.0157;
		TR << "pdbnum test_resnum_enzdes_csts has a result of " << compute_result << ", compared to expected result of " << expect_result << std::endl;
		TS_ASSERT_DELTA( compute_result, expect_result, TOLERATED_ERROR );

		//Create a second, identical EnzScoreFilter object, this time using pose numbering to specify the resnum.
		EnzScoreFilter posenum_test( "107", default_cstid, scorefxn, default_scoretype, default_threshold, default_wholepose, true );

		//Compute the enzdes-style energy for the given posenum, and compare to the value computed for the pdbnum
		core::Real compute_posenum_result = posenum_test.compute(pose);
		TR << "posenum test_resnum_enzdes_csts has a result of " << compute_posenum_result << ", compared to the previously calculated (pdbnum) result of " << compute_result << std::endl;
		TS_ASSERT_EQUALS( compute_result, compute_posenum_result );

		//Create a third EnzScoreFilter object, this time setting resnum to default and specifying the cstid
		EnzScoreFilter cstid_test( default_resnum, "1A", scorefxn, default_scoretype, default_threshold, default_wholepose, true );

		//Compute the enzdes-style energy for the given cstid, and compare to the value computed for the pdbnum specifier
		core::Real compute_cstid_result = cstid_test.compute(pose);
		TR << "cstid test_resnum_enzdes_csts has a result of " << compute_cstid_result << ", compared to the previously calculated (pdbnum) result of " << compute_result << std::endl;
		TS_ASSERT_EQUALS( compute_result, compute_cstid_result );

		TR << "End of test_resnum_enzdes_csts." << std::endl;
	}

	// Unit test for calculation of cst energy using the fully automatic/default mode.
	// Rosetta attempts to automatically identify the ligand.
	void test_auto_enzdes_csts() {
		using namespace protocols::enzdes;
		using namespace core::scoring;

		//Set up a scorefunction with csts
		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->reset();
		scorefxn->set_weight( scoring::atom_pair_constraint, 1.0);
		scorefxn->set_weight( scoring::angle_constraint, 1.0);
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0);

		//Set up variables containing the defaults for each EnzScoreFilter option.
		std::string default_resnum = EnzScoreFilter::default_value_for_resnum();
		std::string default_cstid = EnzScoreFilter::default_value_for_cstid();
		core::scoring::ScoreType default_scoretype = core::scoring::score_type_from_name( EnzScoreFilter::default_value_for_score_type() );
		core::Real default_threshold = EnzScoreFilter::default_value_for_threshold();
		bool default_wholepose = EnzScoreFilter::default_value_for_whole_pose();

		//Add the constraints to the pose
		cst_io->add_constraints_to_pose( pose, scorefxn, true );

		//Create an EnzScoreFilter object.
		//The "true" argument sets is_cstE to true (which overrides the score_type setting).
		EnzScoreFilter auto_test( default_resnum, default_cstid, scorefxn, default_scoretype, default_threshold, default_wholepose, true );

		//Compute the enzdes-style energy
		core::Real compute_result = auto_test.compute(pose);
		core::Real expect_result = 40.0157;
		TR << "test_auto_enzdes_csts has a result of " << compute_result << ", compared to expected result of " << expect_result << std::endl;
		TS_ASSERT_DELTA( compute_result, expect_result, TOLERATED_ERROR );

		TR << "End of test_auto_enzdes_csts." << std::endl;
	}
};


