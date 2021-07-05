// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/dunbrack/DunbrackEnergyTests_peptoid.cxxtest.hh
/// @brief  Test the Dunbrack energy with peptoids, in which derivatives affect phi, psi, and the PRECEDING omega.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Unit headers
#include <core/pack/dunbrack/DunbrackEnergy.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/PartialAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/MinimizationData.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //To make it easier to build poses

// Utility, etc Headers
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh> // AUTO IWYU For ResidueType

static basic::Tracer TR("DunbrackEnergyTests_peptoid");


class DunbrackEnergyTests_peptoid : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-out:levels core.pack.rotamer_set.ContinuousRotamerSet:500 -extra_res_fa core/pack/dunbrack/dummy_type.params" );
	}

	void tearDown(){

	}

	void test_peptoid_dunbrack_scoring(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, nullptr, "");
		builder.add_residue("Append", "601", 2, false, "N", 0, 1, nullptr, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, nullptr, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_psi(2, -63);
		pose.set_omega(2, 180);
		pose.set_omega(1, 179);
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
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score6;
		{
			core::pose::Pose pose6(pose);
			pose6.set_omega(2, -27);
			score6 = (*sfxn)(pose6);
		}
		core::Real score8;
		{
			core::pose::Pose pose8(pose);
			pose8.set_omega(1, -37);
			score8 = (*sfxn)(pose8);
		}



		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score4 << " " << score6 << " " << score8 << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score4-score1) );
		TS_ASSERT_DELTA( score6, score1, 0.00001 );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score8-score1) );
	}

	void test_peptoid_dunbrack_derivs(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, nullptr, "");
		builder.add_residue("Append", "601", 2, false, "N", 0, 1, nullptr, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, nullptr, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_psi(2, -63);
		pose.set_omega(2, 175);
		pose.set_omega(1, 179);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		(*sfxn)(pose); //Ensure that scoring objects are set up.

		core::pack::dunbrack::DunbrackEnergy dun;

		{
			core::scoring::ResSingleMinimizationData mindata;
			core::Real const score1a (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("N"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 1 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1b (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CA"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 2 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1c (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
			core::Real const score1d (dun.eval_residue_dof_derivative( pose.residue(2), mindata, core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 1 ), core::id::PHI ), core::id::TorsionID( 1, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );

			TR.precision( 10 );
			TR << "Derivatives are: " << score1a << " " << score1b << " " << score1c << " " << score1d << std::endl;

			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1a) );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1b) );
			TS_ASSERT_DELTA( score1c, 0.0, 0.00001 );
			TS_ASSERT_LESS_THAN( 0.1, std::abs(score1d) );
		}
	}


	void test_peptoid_IV_reporting() {
		using core::id::PartialAtomID;

		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, nullptr, "");
		builder.add_residue("Append", "601", 2, false, "N", 0, 1, nullptr, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, nullptr, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		core::pack::dunbrack::DunbrackEnergy dunE;
		utility::vector1< PartialAtomID > partial_ids = dunE.atoms_with_dof_derivatives( pose.residue(2), pose );

		TS_ASSERT_EQUALS( partial_ids.size(), 6 );
		TS_ASSERT_EQUALS( partial_ids[1], PartialAtomID( 1, 1, 0 ));
		TS_ASSERT_EQUALS( partial_ids[2], PartialAtomID( 1, 1, 1 ));
		TS_ASSERT_EQUALS( partial_ids[3], PartialAtomID( 1, 2 ));
		TS_ASSERT_EQUALS( partial_ids[4], PartialAtomID( 2, 2 ));
		TS_ASSERT_EQUALS( partial_ids[5], PartialAtomID( 3, 2 ));
		TS_ASSERT_EQUALS( partial_ids[6], PartialAtomID( 1, 3, 0 ));

	}

};
