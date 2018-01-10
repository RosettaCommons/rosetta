// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/dunbrack/DunbrackEnergyTests_oligourea.cxxtest.hh
/// @brief  Test the Dunbrack energy with a residue type in which only some of the mainchain dihedrals affect rotamers (oligoureas).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <numeric/xyzVector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //To make it easier to build poses
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/relax/FastRelax.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("DunbrackEnergyTests_oligourea");


class DunbrackEnergyTests_oligourea : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-out:levels core.pack.rotamer_set.ContinuousRotamerSet:500 -extra_res_fa core/pack/dunbrack/dummy_type.params -output_virtual" );
	}

	void tearDown(){

	}

	void test_beta_aa_dunbrack_scoring(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "B3C", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);
		pose.set_chi(2, 2, 62.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_omega(2, -27);
			score5 = (*sfxn)(pose5);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << std::endl;

		TS_ASSERT_LESS_THAN( 0.01, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score3-score1) );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score4-score1) );
		TS_ASSERT_DELTA( score5, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_scoring(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_VAL", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_mu(2, -77);
			score5 = (*sfxn)(pose5);
		}
		core::Real score6;
		{
			core::pose::Pose pose6(pose);
			pose6.set_omega(2, -27);
			score6 = (*sfxn)(pose6);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << " " << score6 << std::endl;

		TS_ASSERT_LESS_THAN( 0.01, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score3-score1) );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score4-score1) );
		TS_ASSERT_DELTA( score5, score1, 0.00001 );
		TS_ASSERT_DELTA( score6, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_scoring2(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_DUMMYTYPE", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_mu(2, -77);
			score5 = (*sfxn)(pose5);
		}
		core::Real score6;
		{
			core::pose::Pose pose6(pose);
			pose6.set_omega(2, -27);
			score6 = (*sfxn)(pose6);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << " " << score6 << std::endl;

		TS_ASSERT_LESS_THAN( 0.01, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score3-score1) );
		TS_ASSERT_DELTA( score4, score1, 0.00001 );
		TS_ASSERT_LESS_THAN( 0.01, std::abs(score5-score1) );
		TS_ASSERT_DELTA( score6, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_derivs(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_VAL", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, -175);
		pose.set_omega(2, 175);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		(*sfxn)(pose); //Ensure that scoring objects are set up.

		core::pack::dunbrack::DunbrackEnergy dun;

		{
			core::scoring::ResSingleMinimizationData mindata;
			core::Real const score1a (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("N"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 1 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1b (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CA"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 2 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1c (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CM"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1d (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("NU"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 4 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1e (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 5 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );

			TR.precision( 10 );
			TR << "Derivatives are: " << score1a << " " << score1b << " " << score1c << " " << score1d << " " << score1e << std::endl;

			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1a) );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1b) );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1c) );
			TS_ASSERT_DELTA( score1d, 0.0, 0.00001 );
			TS_ASSERT_DELTA( score1e, 0.0, 0.00001 );
		}
	}

	void test_oligourea_dunbrack_derivs2(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_DUMMYTYPE", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, -175);
		pose.set_omega(2, 175);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		(*sfxn)(pose); //Ensure that scoring objects are set up.

		core::pack::dunbrack::DunbrackEnergy dun;

		{
			core::scoring::ResSingleMinimizationData mindata;
			core::Real const score1a (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("N"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 1 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1b (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CA"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 2 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1c (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CM"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1d (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("NU"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 4 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1e (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 5 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );

			TR.precision( 10 );
			TR << "Derivatives are: " << score1a << " " << score1b << " " << score1c << " " << score1d << " " << score1e << std::endl;

			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1a) );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1b) );
			TS_ASSERT_DELTA( score1c, 0.0, 0.00001 );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1d) );
			TS_ASSERT_DELTA( score1e, 0.0, 0.00001 );
		}
	}

	/// @brief Test fa_dun scoring of oligourea-proline.
	void test_oligourea_proline_score() {
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "ALA:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_PRO", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "ALA:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		pose.set_phi(1, -135);
		pose.set_psi(1, 135);
		pose.set_omega(1, 180);
		pose.set_phi(3, -135);
		pose.set_psi(3, 135);
		pose.set_omega(3, 180);

		//Test case:
		pose.set_phi(2, -120);
		pose.set_theta(2, -150);
		pose.set_psi(2, 50);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, -25.1);
		pose.set_chi(2, 2, 36.2);
		pose.set_chi(3, 2, -31.7);
		pose.update_residue_neighbors();
		core::Real const rot1score( (*sfxn)(pose) );
		pose.set_chi(1, 2, 28.5);
		pose.set_chi(2, 2, -34.7);
		pose.set_chi(3, 2, 26.4);
		pose.update_residue_neighbors();
		core::Real const rot2score( (*sfxn)(pose) );
		pose.set_chi(1, 2, 77);
		pose.set_chi(2, 2, -94.7);
		pose.set_chi(3, 2, 102.4);
		pose.update_residue_neighbors();
		core::Real const nonrotscore( (*sfxn)(pose) );

		TS_ASSERT_LESS_THAN( rot1score, nonrotscore );
		TS_ASSERT_LESS_THAN( rot2score, nonrotscore );
		TS_ASSERT_LESS_THAN( rot1score, rot2score );
		TS_ASSERT_LESS_THAN( rot1score, 10.0 );
		TS_ASSERT_LESS_THAN( rot2score, 10.0 );

		TR << "Rot1Score:\t" << rot1score << std::endl;
		TR << "Rot2Score:\t" << rot2score << std::endl;
		TR << "NonRotScore:\t" << nonrotscore << std::endl;
	}

	/// @brief Test FastRelax with oligourea-proline, confirming that pro_close keeps the proline ring closed.
	void test_oligourea_proline_pack() {
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "ALA:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_PRO", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "ALA:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(1, -135);
		pose.set_psi(1, 135);
		pose.set_omega(1, 180);
		pose.set_phi(2, -77.6);
		pose.set_theta(2, 48.1);
		pose.set_psi(2, 102.2);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_phi(3, -135);
		pose.set_psi(3, 135);
		pose.set_omega(3, 180);
		pose.set_chi(1, 2, 100); //Set a ridiculous chi
		pose.set_chi(2, 2, 100); //Set a ridiculous chi
		pose.set_chi(3, 2, 100); //Set a ridiculous chi
		pose.update_residue_neighbors();

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );
		sfxn->set_weight( core::scoring::pro_close, 1.0 );

		(*sfxn)(pose);

		core::Real const pro_dist1( pose.xyz( core::id::NamedAtomID("NV", 2) ).distance( pose.xyz( core::id::NamedAtomID("N", 2) ) ) );
		TR << "Pre-relax proline distance: " << pro_dist1 << std::endl;

		TS_ASSERT_LESS_THAN( 1.0, pro_dist1 );

		protocols::relax::FastRelax frlx(sfxn, 1);
		frlx.apply(pose);

		/*core::pack::task::TaskFactoryOP taskfact( new core::pack::task::TaskFactory );
		taskfact->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::RestrictToRepacking ) );
		protocols::simple_moves::PackRotamersMover pack(sfxn, taskfact->create_task_and_apply_taskoperations(pose));
		pack.apply(pose);*/

		core::Real const pro_dist2( pose.xyz( core::id::NamedAtomID("NV", 2) ).distance( pose.xyz( core::id::NamedAtomID("N", 2) ) ) );
		TR << "Post-relax proline distance: " << pro_dist2 << std::endl;

		TS_ASSERT_LESS_THAN( pro_dist2, pro_dist1 );
		TS_ASSERT_LESS_THAN( pro_dist2, 0.1 );

	}

};
