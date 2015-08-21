// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/FavorSequenceProfile.cxxtest.hh
/// @brief  test for FavorSequeceProfile mover
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_moves/FavorSequenceProfile.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.simple_moves.FavorSequenceProfile.cxxtest.hh");

// --------------- Test Class --------------- //

class FavorSequenceProfileTests : public CxxTest::TestSuite {

private:
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();

		scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		scorefxn_->set_weight(core::scoring::res_type_constraint, 1.0);
	}

	void tearDown() {
	}

	void test_use_IDENTITY() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		core::sequence::Sequence seq(pose);
		FavorSequenceProfile FSP;
		FSP.set_sequence(seq,"IDENTITY");
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have an energy bonus equal to weight * (number of residues)
		core::Real favored_score = scorefxn_->score(pose);
		TS_ASSERT_DELTA(favored_score, (-1.0 * pose.total_residue()), 0.01 );

		//Spotcheck that mis-match causes zero score
		TS_ASSERT_EQUALS( pose.aa(2), core::chemical::aa_ala );
		core::conformation::Residue res_ala(pose.residue_type(2), true); // The true is a dummy arg.
		TS_ASSERT_EQUALS( pose.aa(1), core::chemical::aa_asp ); // 49-48 = 1
		pose.replace_residue(1, res_ala, true);
		TS_ASSERT_EQUALS( pose.aa(11), core::chemical::aa_trp ); // 59-48 = 11
		pose.replace_residue(11, res_ala, true);
		TS_ASSERT_EQUALS( pose.aa(60), core::chemical::aa_ile ); //108-48 = 60
		pose.replace_residue(60, res_ala, true);
		TS_ASSERT_EQUALS( pose.aa(100), core::chemical::aa_lys ); // 148-48 = 100
		pose.replace_residue(100, res_ala, true);

		scorefxn_->score(pose);
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(1)[core::scoring::res_type_constraint], 0, 0.0001);
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(11)[core::scoring::res_type_constraint], 0, 0.0001);
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(60)[core::scoring::res_type_constraint], 0, 0.0001);
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(100)[core::scoring::res_type_constraint], 0, 0.0001);
	}

	void test_use_weight() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		core::sequence::Sequence seq(pose);
		FavorSequenceProfile FSP;
		FSP.set_weight( 0.35 );
		FSP.set_sequence(seq,"IDENTITY");
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have an energy bonus equal to weight * (number of residues)
		core::Real favored_score = scorefxn_->score(pose);
		TS_ASSERT_DELTA( favored_score , (-0.35 * pose.total_residue()) , 0.01);
	}

	void test_use_BLOSUM62() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		core::sequence::Sequence seq(pose);
		FavorSequenceProfile FSP;
		FSP.set_sequence(seq,"BLOSUM62");
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have a net energy bonus.
		core::Real favored_score = scorefxn_->score(pose);
		TR << "Blosum62 total: " << favored_score << std::endl;
		TS_ASSERT_LESS_THAN(favored_score, 0);

		// Check each individual residues value
		for ( core::Size ii(1); ii <= pose.total_residue(); ++ii ) {
			core::Real resenergy = pose.energies().residue_total_energies(ii)[core::scoring::res_type_constraint];
			//TR << "Blosum62 resi " << ii << " res_type_constraint score " << resenergy << std::endl;
			TS_ASSERT_LESS_THAN(resenergy, 0);
		}

		// Spot check individual positions
		// The exact values here aren't important, merely that they match up with what they should be from the matrix (aa/pos wise)
		TS_ASSERT_EQUALS( pose.aa(1), core::chemical::aa_asp ); // 49-48 = 1
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(1)[core::scoring::res_type_constraint], -0.9652, 0.0001); //Values are boltzman adjusted
		TS_ASSERT_EQUALS( pose.aa(11), core::chemical::aa_trp ); // 59-48 = 11
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(11)[core::scoring::res_type_constraint], -0.9998, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(60), core::chemical::aa_ile ); //108-48 = 60
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(60)[core::scoring::res_type_constraint], -0.6214, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(100), core::chemical::aa_lys ); // 148-48 = 100
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(100)[core::scoring::res_type_constraint], -0.8931, 0.0001);

	}

	void test_use_PSSM() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		std::string pssm_name("protocols/simple_moves/1bl0_auto.fasta.pssm");
		core::sequence::SequenceProfile profile;
		profile.read_from_file(pssm_name);
		FavorSequenceProfile FSP;
		FSP.set_profile(profile);
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have a net energy bonus.
		core::Real favored_score = scorefxn_->score(pose);
		TR << "PSSM total: " << favored_score << std::endl;
		TS_ASSERT_LESS_THAN(favored_score, 0);

		// Check each individual residues value
		for ( core::Size ii(1); ii <= pose.total_residue(); ++ii ) {
			core::Real resenergy = pose.energies().residue_total_energies(ii)[core::scoring::res_type_constraint];
			//TR << "Blosum62 resi " << ii << " res_type_constraint score " << resenergy << std::endl;
			TS_ASSERT_LESS_THAN(resenergy, 0);
		}

		// Spot check individual positions
		// The exact values here aren't important, merely that they match up with what they should be from the matrix (aa/pos wise)
		TS_ASSERT_EQUALS( pose.aa(1), core::chemical::aa_asp ); // 49-48 = 1
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(1)[core::scoring::res_type_constraint], -0.9790, 0.0001); //Values are boltzman adjusted
		TS_ASSERT_EQUALS( pose.aa(11), core::chemical::aa_trp ); // 59-48 = 11
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(11)[core::scoring::res_type_constraint], -0.4976, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(60), core::chemical::aa_ile ); //108-48 = 60
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(60)[core::scoring::res_type_constraint], -0.0009, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(100), core::chemical::aa_lys ); // 148-48 = 100
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(100)[core::scoring::res_type_constraint], -0.0309, 0.0001);

	}

	void test_no_scaling() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		core::sequence::Sequence seq(pose);
		FavorSequenceProfile FSP;
		FSP.set_sequence(seq,"BLOSUM62");
		FSP.set_scaling( "none" );
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have a net energy bonus.
		core::Real favored_score = scorefxn_->score(pose);
		TR << "Blosum62 total: " << favored_score << std::endl;
		TS_ASSERT_LESS_THAN(favored_score, 0);

		// Check each individual residues value
		for ( core::Size ii(1); ii <= pose.total_residue(); ++ii ) {
			core::Real resenergy = pose.energies().residue_total_energies(ii)[core::scoring::res_type_constraint];
			//TR << "Blosum62 resi " << ii << " res_type_constraint score " << resenergy << std::endl;
			TS_ASSERT_LESS_THAN(resenergy, 0);
		}

		// Spot check individual positions
		// The exact values here aren't important, merely that they match up with what they should be from the matrix (aa/pos wise)
		TS_ASSERT_EQUALS( pose.aa(1), core::chemical::aa_asp ); // 49-48 = 1
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(1)[core::scoring::res_type_constraint], -6.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(11), core::chemical::aa_trp ); // 59-48 = 11
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(11)[core::scoring::res_type_constraint], -11.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(60), core::chemical::aa_ile ); //108-48 = 60
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(60)[core::scoring::res_type_constraint], -4.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(100), core::chemical::aa_lys ); // 148-48 = 100
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(100)[core::scoring::res_type_constraint], -5.0, 0.0001);
	}

	void test_global_scaling() {
		using namespace protocols::simple_moves;

		core::pose::Pose pose(create_test_in_pdb_pose());
		core::Real original_score = scorefxn_->score(pose);

		core::sequence::Sequence seq(pose);
		FavorSequenceProfile FSP;
		FSP.set_sequence(seq,"BLOSUM62");
		FSP.set_scaling( "global" );
		FSP.apply(pose);

		TS_ASSERT_EQUALS(original_score, 0);
		// All residues match, so we should have a net energy bonus.
		core::Real favored_score = scorefxn_->score(pose);
		TR << "Blosum62 total: " << favored_score << std::endl;
		TS_ASSERT_LESS_THAN(favored_score, 0);

		// Check each individual residues value
		for ( core::Size ii(1); ii <= pose.total_residue(); ++ii ) {
			core::Real resenergy = pose.energies().residue_total_energies(ii)[core::scoring::res_type_constraint];
			//TR << "Blosum62 resi " << ii << " res_type_constraint score " << resenergy << std::endl;
			TS_ASSERT_LESS_THAN(resenergy, 0);
		}

		// Spot check individual positions
		// The exact values here aren't important, merely that they match up with what they should be from the matrix (aa/pos wise)
		TS_ASSERT_EQUALS( pose.aa(1), core::chemical::aa_asp ); // 49-48 = 1
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(1)[core::scoring::res_type_constraint], -6.0/11.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(11), core::chemical::aa_trp ); // 59-48 = 11
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(11)[core::scoring::res_type_constraint], -11.0/11.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(60), core::chemical::aa_ile ); //108-48 = 60
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(60)[core::scoring::res_type_constraint], -4.0/11.0, 0.0001);
		TS_ASSERT_EQUALS( pose.aa(100), core::chemical::aa_lys ); // 148-48 = 100
		TS_ASSERT_DELTA( pose.energies().residue_total_energies(100)[core::scoring::res_type_constraint], -5.0/11.0, 0.0001);

	}


};
