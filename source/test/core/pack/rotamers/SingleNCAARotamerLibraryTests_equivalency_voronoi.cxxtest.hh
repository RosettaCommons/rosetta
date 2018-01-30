// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/rotamers/SingleNCAARotamerLibraryTests_equivalency_voronoi.cxxtest.hh
/// @brief  Unit tests for the non-canonical rotamer libraries.  These tests ensure that canonical and non-canonical fa_dun
/// scoring is identical, with Voronoi-based detection of rotamer wells.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleNCAARotamerLibraryCreator.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SingleNCAARotamerLibraryTests_equivalency_voronoi");


class SingleNCAARotamerLibraryTests_equivalency_voronoi : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		// This is just leucine, but treated as though it's a non-canonical amino acid.
		core_init_with_additional_options("-extra_res_fa core/pack/rotamers/LEU_NCAA.params core/pack/rotamers/DAB_corrected_order.params -fa_dun_canonicals_use_voronoi true -fa_dun_noncanonicals_use_voronoi true");
	}

	void tearDown(){

	}

	/// @author Andy Watkins
	void test_equivalent_rotlib_values_to_canonical() {
		using namespace core;
		using namespace core::chemical;
		using namespace core::pack::dunbrack;
		using namespace core::pack::rotamers;

		pose::Pose pose_ca, pose_ncaa;
		pose::make_pose_from_sequence(pose_ca, "GLG", "fa_standard");
		pose::make_pose_from_sequence(pose_ncaa, "GX[LEU_NCAA]G", "fa_standard");

		auto rlf = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance();
		using namespace core::pack::dunbrack;

		SingleResidueDunbrackLibraryCOP rl_cleu = utility::pointer::dynamic_pointer_cast< RotamericSingleResidueDunbrackLibrary< TWO, TWO > const >( rlf->get( pose_ca.residue_type( 2 ) ) );
		SingleResidueDunbrackLibraryCOP rl_nleu = utility::pointer::dynamic_pointer_cast< RotamericSingleResidueDunbrackLibrary< TWO, TWO > const >( rlf->get( pose_ncaa.residue_type( 2 ) ) );

		for ( Real phi = -180; phi <= 180; phi += 5 ) {
			for ( Real psi = -180; psi <= 180; psi += 5 ) {
				for ( Size ii = 1; ii <= 9; ++ii ) {
					TS_ASSERT( rl_cleu->get_probability_for_rotamer( phi, psi, ii ) == rl_nleu->get_probability_for_rotamer( phi, psi, ii ) );
				}
			}
		}


		// Are we getting the same rotamer samples, in the same order?
		for ( Real phi = -180; phi <= 180; phi += 1 ) {
			for ( Real psi = -180; psi <= 180; psi += 1 ) {
				utility::vector1< DunbrackRotamerSampleData > nleu_samples = rl_nleu->get_all_rotamer_samples( utility::fixedsizearray1< Real, 5 >({phi, psi,0,0,0}) );
				utility::vector1< DunbrackRotamerSampleData > cleu_samples = rl_cleu->get_all_rotamer_samples( utility::fixedsizearray1< Real, 5 >({phi, psi,0,0,0}) );
				TS_ASSERT( nleu_samples.size() == cleu_samples.size() );
				for ( Size ii = 1; ii <= cleu_samples.size(); ++ii ) {
					TS_ASSERT( nleu_samples[ii].nrchi_sample() == cleu_samples[ii].nrchi_sample() );
					TS_ASSERT( nleu_samples[ii].nchi() == cleu_samples[ii].nchi() );
					TS_ASSERT( nleu_samples[ii].rot_well() == cleu_samples[ii].rot_well() );
					TS_ASSERT( nleu_samples[ii].chi_mean() == cleu_samples[ii].chi_mean() );
					TS_ASSERT( nleu_samples[ii].chi_sd() == cleu_samples[ii].chi_sd() );
					TS_ASSERT( nleu_samples[ii].probability() == cleu_samples[ii].probability() );
				}
			}
		}

		// What about energies?
		for ( Real phi = -180; phi <= 180; phi += 15 ) {
			for ( Real psi = -180; psi <= 180; psi += 15 ) {
				for ( Real chi1 = -10; chi1 <= 130; chi1 += 15 ) {
					for ( Real chi2 = -10; chi2 <= 130; chi2 += 15 ) {
						RotamerLibraryScratchSpace nscratch, cscratch;
						pose_ca.set_phi( 2, phi );
						pose_ncaa.set_phi( 2, phi );
						pose_ca.set_psi( 2, psi );
						pose_ncaa.set_psi( 2, psi );
						pose_ca.set_chi( 1, 2, chi1 );
						pose_ca.set_chi( 2, 2, chi2 );
						pose_ncaa.set_chi( 1, 2, chi1 );
						pose_ncaa.set_chi( 2, 2, chi2 );
						Real cv = rl_cleu->rotamer_energy( pose_ca.residue( 2 ), pose_ca, cscratch );
						Real nv = rl_nleu->rotamer_energy( pose_ncaa.residue( 2 ), pose_ncaa, nscratch );
						TS_ASSERT( cv == nv );
						TS_ASSERT( cscratch.chimean() == nscratch.chimean() );
						TS_ASSERT( cscratch.chisd() == nscratch.chisd() );
						TS_ASSERT( cscratch.chidev() == nscratch.chidev() );
						TS_ASSERT( cscratch.chidevpen() == nscratch.chidevpen() );
						TS_ASSERT( cscratch.dchidevpen_dchi() == nscratch.dchidevpen_dchi() );
						TS_ASSERT( cscratch.dchidevpen_dbb() == nscratch.dchidevpen_dbb() );
						TS_ASSERT( cscratch.drotprob_dbb() == nscratch.drotprob_dbb() );
						TS_ASSERT( cscratch.dneglnrotprob_dbb() == nscratch.dneglnrotprob_dbb() );
						TS_ASSERT( cscratch.dchimean_dbb() == nscratch.dchimean_dbb() );
						TS_ASSERT( cscratch.dchisd_dbb() == nscratch.dchisd_dbb() );

					}
				}
			}
		}
	}

	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_equivalent_scoring_to_canonical() {
		core::pose::Pose pose_ca, pose_ncaa;
		core::pose::make_pose_from_sequence(pose_ca, "GLG", "fa_standard");
		core::pose::make_pose_from_sequence(pose_ncaa, "GX[LEU_NCAA]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(5); //1 = phi, 2=psi, 3=chi1, 4=chi2
		trial_conformations[1][1] = -61.0; trial_conformations[1][2] = -41.0; trial_conformations[1][3] = 32.0; trial_conformations[1][4] = -53.0;
		trial_conformations[2][1] = -135.0; trial_conformations[2][2] = 135.0; trial_conformations[2][3] = 32.0; trial_conformations[2][4] = -53.0;
		trial_conformations[3][1] = -61.0; trial_conformations[3][2] = -41.0; trial_conformations[3][3] = 72.0; trial_conformations[3][4] = 63.0;
		trial_conformations[4][1] = -61.0; trial_conformations[4][2] = -41.0; trial_conformations[4][3] = -42.0; trial_conformations[4][4] = -122.0;
		trial_conformations[5][1] = -135.0; trial_conformations[5][2] = 135.0; trial_conformations[5][3] = -42.0; trial_conformations[5][4] = -122.0;

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tPSI\tCHI1\tCHI2\tSCORE_CA\tSCORE_NCAA" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose_ca.set_phi( 2, trial_conformations[i][1] );
			pose_ca.set_psi( 2, trial_conformations[i][2] );
			pose_ca.set_chi( 1, 2, trial_conformations[i][3] );
			pose_ca.set_chi( 2, 2, trial_conformations[i][4] );
			pose_ncaa.set_phi( 2, trial_conformations[i][1] );
			pose_ncaa.set_psi( 2, trial_conformations[i][2] );
			pose_ncaa.set_chi( 1, 2, trial_conformations[i][3] );
			pose_ncaa.set_chi( 2, 2, trial_conformations[i][4] );
			core::Real const score_ca( sfxn(pose_ca) );
			core::Real const score_ncaa( sfxn(pose_ncaa) );
			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << score_ca << "\t" << score_ncaa << std::endl;
			TS_ASSERT_DELTA( score_ca, score_ncaa, 1e-6 );
		}

		TR << std::endl;
	}

	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_equivalent_derivatives_to_canonical() {
		core::pose::Pose pose_ca, pose_ncaa;
		core::pose::make_pose_from_sequence(pose_ca, "GLG", "fa_standard");
		core::pose::make_pose_from_sequence(pose_ncaa, "GX[LEU_NCAA]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(5); //1 = phi, 2=psi, 3=chi1, 4=chi2
		trial_conformations[1][1] = -61.0; trial_conformations[1][2] = -41.0; trial_conformations[1][3] = 32.0; trial_conformations[1][4] = -53.0;
		trial_conformations[2][1] = -135.0; trial_conformations[2][2] = 135.0; trial_conformations[2][3] = 32.0; trial_conformations[2][4] = -53.0;
		trial_conformations[3][1] = -61.0; trial_conformations[3][2] = -41.0; trial_conformations[3][3] = 72.0; trial_conformations[3][4] = 63.0;
		trial_conformations[4][1] = -61.0; trial_conformations[4][2] = -41.0; trial_conformations[4][3] = -42.0; trial_conformations[4][4] = -122.0;
		trial_conformations[5][1] = -135.0; trial_conformations[5][2] = 135.0; trial_conformations[5][3] = -42.0; trial_conformations[5][4] = -122.0;

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		core::pack::dunbrack::DunbrackEnergy dun;

		TR << "PHI\tPSI\tCHI1\tCHI2\tDERIV_CA_PHI\tDERIV_NCAA_PHI\tDERIV_CA_PSI\tDERIV_NCAA_PSI\tDERIV_CA_CHI1\tDERIV_NCAA_CHI1\tDERIV_CA_CHI2\tDERIV_NCAA_CHI2\tDERIV_CA_OMG\tDERIV_NCAA_OMG" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose_ca.set_phi( 2, trial_conformations[i][1] );
			pose_ca.set_psi( 2, trial_conformations[i][2] );
			pose_ca.set_chi( 1, 2, trial_conformations[i][3] );
			pose_ca.set_chi( 2, 2, trial_conformations[i][4] );
			pose_ncaa.set_phi( 2, trial_conformations[i][1] );
			pose_ncaa.set_psi( 2, trial_conformations[i][2] );
			pose_ncaa.set_chi( 1, 2, trial_conformations[i][3] );
			pose_ncaa.set_chi( 2, 2, trial_conformations[i][4] );

			core::scoring::ResSingleMinimizationData rmm;
			core::Real const dEdphi_ca(    dun.eval_residue_dof_derivative( pose_ca.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 1), pose_ca, sfxn, sfxn.weights() ) );
			core::Real const dEdphi_ncaa(  dun.eval_residue_dof_derivative( pose_ncaa.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 1), pose_ncaa, sfxn, sfxn.weights() ) );
			core::Real const dEdpsi_ca(    dun.eval_residue_dof_derivative( pose_ca.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 2), pose_ca, sfxn, sfxn.weights() ) );
			core::Real const dEdpsi_ncaa(  dun.eval_residue_dof_derivative( pose_ncaa.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 2), pose_ncaa, sfxn, sfxn.weights() ) );
			core::Real const dEdchi1_ca(   dun.eval_residue_dof_derivative( pose_ca.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::CHI, 1), pose_ca, sfxn, sfxn.weights() ) );
			core::Real const dEdchi1_ncaa( dun.eval_residue_dof_derivative( pose_ncaa.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::CHI, 1), pose_ncaa, sfxn, sfxn.weights() ) );
			core::Real const dEdchi2_ca(   dun.eval_residue_dof_derivative( pose_ca.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::CHI, 2), pose_ca, sfxn, sfxn.weights() ) );
			core::Real const dEdchi2_ncaa( dun.eval_residue_dof_derivative( pose_ncaa.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::CHI, 2), pose_ncaa, sfxn, sfxn.weights() ) );
			core::Real const dEdomg_ca(    dun.eval_residue_dof_derivative( pose_ca.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 3), pose_ca, sfxn, sfxn.weights() ) );
			core::Real const dEdomg_ncaa(  dun.eval_residue_dof_derivative( pose_ncaa.residue( 2 ), rmm, core::id::DOF_ID( core::id::AtomID( 1, 2 ), core::id::PHI ), core::id::TorsionID(2, core::id::BB, 3), pose_ncaa, sfxn, sfxn.weights() ) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << dEdphi_ca << "\t" << dEdphi_ncaa << "\t" << dEdpsi_ca << "\t" << dEdpsi_ncaa << "\t" << dEdchi1_ca << "\t" << dEdchi1_ncaa << "\t" << dEdchi2_ca << "\t" << dEdchi2_ncaa << "\t" << dEdomg_ca << "\t" << dEdomg_ncaa << std::endl;
			TS_ASSERT_DELTA( dEdphi_ca, dEdphi_ncaa, 1e-6 );
			TS_ASSERT_DELTA( dEdpsi_ca, dEdpsi_ncaa, 1e-6 );
			TS_ASSERT_DELTA( dEdchi1_ca, dEdchi1_ncaa, 1e-6 );
			TS_ASSERT_DELTA( dEdchi2_ca, dEdchi2_ncaa, 1e-6 );
			TS_ASSERT_DELTA( dEdomg_ca, dEdomg_ncaa, 1e-6 );
			TS_ASSERT_DELTA( dEdomg_ca, 0.0, 1e-6 );
			TS_ASSERT_DELTA( dEdomg_ncaa, 0.0, 1e-6 );
		}

		TR << std::endl;
	}

};
