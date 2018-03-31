// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/rotamers/SingleNCAARotamerLibraryTests.cxxtest.hh
/// @brief  Unit tests for the non-canonical rotamer libraries.
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

static basic::Tracer TR("SingleNCAARotamerLibraryTests");


class SingleNCAARotamerLibraryTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		// This is just leucine, but treated as though it's a non-canonical amino acid.
		core_init_with_additional_options("-extra_res_fa core/pack/rotamers/LEU_NCAA.params core/pack/rotamers/DAB_corrected_order.params");
	}

	void tearDown(){

	}

	/// @brief Confirm that the rotamer wells in the B53 (DAB) rotamer library really do have
	/// the appropriate score values.
	void test_DAB_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[DAB]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(9); //1 = phi, 2=psi, 3=chi1, 4=chi2
		trial_conformations[1][1] = -60.0; trial_conformations[1][2] = -40.0; trial_conformations[1][3] =  -64.0; trial_conformations[1][4] =  178.8;
		trial_conformations[2][1] = -60.0; trial_conformations[2][2] = -40.0; trial_conformations[2][3] =  -62.6; trial_conformations[2][4] =  -64.7;
		trial_conformations[3][1] = -60.0; trial_conformations[3][2] = -40.0; trial_conformations[3][3] = -176.2; trial_conformations[3][4] =  175.8;
		trial_conformations[4][1] = -60.0; trial_conformations[4][2] = -40.0; trial_conformations[4][3] = -176.6; trial_conformations[4][4] =   62.7;
		trial_conformations[5][1] = -60.0; trial_conformations[5][2] = -40.0; trial_conformations[5][3] =  -69.6; trial_conformations[5][4] =   77.5;
		trial_conformations[6][1] = -60.0; trial_conformations[6][2] = -40.0; trial_conformations[6][3] =   40.0; trial_conformations[6][4] = -179.0;
		trial_conformations[7][1] = -60.0; trial_conformations[7][2] = -40.0; trial_conformations[7][3] = -170.7; trial_conformations[7][4] =  -89.1;
		trial_conformations[8][1] = -60.0; trial_conformations[8][2] = -40.0; trial_conformations[8][3] =   34.9; trial_conformations[8][4] =  -85.9;
		trial_conformations[9][1] = -60.0; trial_conformations[9][2] = -40.0; trial_conformations[9][3] =   44.2; trial_conformations[9][4] =  102.4;

		utility::vector1< core::Real > expected_energies( 9 );
		expected_energies[1] = -std::log(0.401037);
		expected_energies[2] = -std::log(0.299396);
		expected_energies[3] = -std::log(0.171242);
		expected_energies[4] = -std::log(0.098707);
		expected_energies[5] = -std::log(0.022515);
		expected_energies[6] = -std::log(0.005180);
		expected_energies[7] = -std::log(0.001878);
		expected_energies[8] = -std::log(0.000042);
		expected_energies[9] = -std::log(0.000004);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tPSI\tCHI1\tCHI2\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_psi( 2, trial_conformations[i][2] );
			pose.set_chi( 1, 2, trial_conformations[i][3] );
			pose.set_chi( 2, 2, trial_conformations[i][4] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

	/// @brief Confirm that the rotamer wells in the ornithine rotamer library really do have
	/// the appropriate score values.
	void test_ORN_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[ORN]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 5 > > trial_conformations(27); //1 = phi, 2=psi, 3=chi1, 4=chi2
		for ( core::Size i(1); i<=trial_conformations.size(); ++i ) {
			trial_conformations[i][1] = -60;
			trial_conformations[i][2] = -40;
		}

		trial_conformations[ 1][3] =  -62.4; trial_conformations[ 1][4] =  179.7; trial_conformations[ 1][5] = -179.8;
		trial_conformations[ 2][3] =  -61.1; trial_conformations[ 2][4] = -178.5; trial_conformations[ 2][5] =  -65.2;
		trial_conformations[ 3][3] =  -63.1; trial_conformations[ 3][4] =  178.1; trial_conformations[ 3][5] =   66.6;
		trial_conformations[ 4][3] =  -61.6; trial_conformations[ 4][4] =  -67.4; trial_conformations[ 4][5] = -178.3;
		trial_conformations[ 5][3] = -174.5; trial_conformations[ 5][4] =  177.5; trial_conformations[ 5][5] = -179.1;
		trial_conformations[ 6][3] = -175.3; trial_conformations[ 6][4] =  175.0; trial_conformations[ 6][5] =   65.3;
		trial_conformations[ 7][3] = -175.5; trial_conformations[ 7][4] =  179.4; trial_conformations[ 7][5] =  -66.6;
		trial_conformations[ 8][3] =  -62.4; trial_conformations[ 8][4] =  -67.9; trial_conformations[ 8][5] =  -63.7;
		trial_conformations[ 9][3] = -174.9; trial_conformations[ 9][4] =   67.9; trial_conformations[ 9][5] =  177.5;
		trial_conformations[10][3] = -174.7; trial_conformations[10][4] =   67.0; trial_conformations[10][5] =   64.3;
		trial_conformations[11][3] =  -72.8; trial_conformations[11][4] =   84.2; trial_conformations[11][5] =  177.9;
		trial_conformations[12][3] =   38.0; trial_conformations[12][4] =  179.7; trial_conformations[12][5] = -179.8;
		trial_conformations[13][3] =  -65.9; trial_conformations[13][4] =   92.3; trial_conformations[13][5] =   65.6;
		trial_conformations[14][3] =   38.1; trial_conformations[14][4] = -178.5; trial_conformations[14][5] =  -65.4;
		trial_conformations[15][3] =   36.7; trial_conformations[15][4] =  177.5; trial_conformations[15][5] =   65.5;
		trial_conformations[16][3] =  -60.7; trial_conformations[16][4] =  -70.4; trial_conformations[16][5] =   90.4;
		trial_conformations[17][3] = -175.3; trial_conformations[17][4] = -105.0; trial_conformations[17][5] = -180.0;
		trial_conformations[18][3] =  -64.2; trial_conformations[18][4] =  121.6; trial_conformations[18][5] =  -64.0;
		trial_conformations[19][3] = -174.2; trial_conformations[19][4] = -106.9; trial_conformations[19][5] =  -66.0;
		trial_conformations[20][3] = -170.8; trial_conformations[20][4] =   97.2; trial_conformations[20][5] =  -69.7;
		trial_conformations[21][3] =   34.8; trial_conformations[21][4] =  -94.9; trial_conformations[21][5] =  177.5;
		trial_conformations[22][3] =   39.8; trial_conformations[22][4] =  135.2; trial_conformations[22][5] = -179.8;
		trial_conformations[23][3] =   40.0; trial_conformations[23][4] =  135.4; trial_conformations[23][5] =   65.2;
		trial_conformations[24][3] =   33.8; trial_conformations[24][4] =  -94.4; trial_conformations[24][5] =  -69.0;
		trial_conformations[25][3] =   42.7; trial_conformations[25][4] = -120.5; trial_conformations[25][5] =   67.6;
		trial_conformations[26][3] = -106.1; trial_conformations[26][4] =  -82.8; trial_conformations[26][5] =   82.9;
		trial_conformations[27][3] =   48.6; trial_conformations[27][4] =  103.1; trial_conformations[27][5] = -151.2;

		utility::vector1< core::Real > expected_energies( 27 );
		expected_energies[ 1] = -std::log(0.248333);
		expected_energies[ 2] = -std::log(0.162028);
		expected_energies[ 3] = -std::log(0.146770);
		expected_energies[ 4] = -std::log(0.115280);
		expected_energies[ 5] = -std::log(0.101497);
		expected_energies[ 6] = -std::log(0.075084);
		expected_energies[ 7] = -std::log(0.053785);
		expected_energies[ 8] = -std::log(0.039650);
		expected_energies[ 9] = -std::log(0.032093);
		expected_energies[10] = -std::log(0.011275);
		expected_energies[11] = -std::log(0.004443);
		expected_energies[12] = -std::log(0.002527);
		expected_energies[13] = -std::log(0.001617);
		expected_energies[14] = -std::log(0.001513);
		expected_energies[15] = -std::log(0.001428);
		expected_energies[16] = -std::log(0.001273);
		expected_energies[17] = -std::log(0.000513);
		expected_energies[18] = -std::log(0.000413);
		expected_energies[19] = -std::log(0.000230);
		expected_energies[20] = -std::log(0.000208);
		expected_energies[21] = -std::log(0.000014);
		expected_energies[22] = -std::log(0.000010);
		expected_energies[23] = -std::log(0.000007);
		expected_energies[24] = -std::log(0.000005);
		expected_energies[25] = -std::log(0.000002);
		expected_energies[26] = -std::log(0.000001);
		expected_energies[27] = -std::log(1e-6);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tPSI\tCHI1\tCHI2\tCHI3\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_psi( 2, trial_conformations[i][2] );
			pose.set_chi( 1, 2, trial_conformations[i][3] );
			pose.set_chi( 2, 2, trial_conformations[i][4] );
			pose.set_chi( 3, 2, trial_conformations[i][5] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << trial_conformations[i][5] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

	/// @brief Confirm that the rotamer wells in the DAP rotamer library really do have
	/// the appropriate score values.
	void test_DPP_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[DPP]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 3 > > trial_conformations(6); //1 = phi, 2=psi, 3=chi1
		for ( core::Size i(1); i<=3; ++i ) {
			trial_conformations[i][1] =  -60;
			trial_conformations[i][2] =  -40;
		}
		for ( core::Size i(4); i<=6; ++i ) {
			trial_conformations[i][1] = -130;
			trial_conformations[i][2] =  130;
		}

		trial_conformations[ 1][3] =  -61.8;
		trial_conformations[ 2][3] = -169.6;
		trial_conformations[ 3][3] =   48.0;
		trial_conformations[ 4][3] =  -63.9;
		trial_conformations[ 5][3] = -169.5;
		trial_conformations[ 6][3] =   50.4;

		utility::vector1< core::Real > expected_energies( 6 );
		expected_energies[ 1] = -std::log(0.838728);
		expected_energies[ 2] = -std::log(0.149772);
		expected_energies[ 3] = -std::log(0.011499);
		expected_energies[ 4] = -std::log(0.788384);
		expected_energies[ 5] = -std::log(0.188797);
		expected_energies[ 6] = -std::log(0.022819);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tPSI\tCHI1\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_psi( 2, trial_conformations[i][2] );
			pose.set_chi( 1, 2, trial_conformations[i][3] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

	/// @brief Confirm that the rotamer wells in the peptoid 601 rotamer library really do have
	/// the appropriate score values.
	void test_peptoid_601_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[601]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 5 > > trial_conformations(6); //1=omega(n-1), 2=phi, 3=psi, 4=chi1, 5=chi2
		for ( core::Size i(1); i<=6; ++i ) {
			trial_conformations[i][1] = -170;
			trial_conformations[i][2] =  -60;
			trial_conformations[i][3] =  -40;
		}

		trial_conformations[ 1][4] =   27.7;
		trial_conformations[ 2][4] =   79.4;
		trial_conformations[ 3][4] =   70.7;
		trial_conformations[ 4][4] = -139.4;
		trial_conformations[ 5][4] = -136.7;
		trial_conformations[ 6][4] =   -1.9;

		trial_conformations[ 1][5] =   74.0;
		trial_conformations[ 2][5] =  102.0;
		trial_conformations[ 3][5] =   28.8;
		trial_conformations[ 4][5] =   49.9;
		trial_conformations[ 5][5] =   41.7;
		trial_conformations[ 6][5] =  -44.0;

		utility::vector1< core::Real > expected_energies( 6 );
		expected_energies[ 1] = -std::log(0.993782);
		expected_energies[ 2] = -std::log(0.006120);
		expected_energies[ 3] = -std::log(0.000098);
		expected_energies[ 4] = -std::log(1e-6);
		expected_energies[ 5] = -std::log(1e-6);
		expected_energies[ 6] = -std::log(1e-6);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "OMEGA\tPHI\tPSI\tCHI1\tCHI2\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_omega( 1, trial_conformations[i][1] );
			pose.set_phi( 2, trial_conformations[i][2] );
			pose.set_psi( 2, trial_conformations[i][3] );
			pose.set_chi( 1, 2, trial_conformations[i][4] );
			pose.set_chi( 2, 2, trial_conformations[i][5] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << trial_conformations[i][5] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

	/// @brief Confirm that the rotamer wells in the B3C rotamer library really do have
	/// the appropriate score values.
	void test_B3C_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[B3C]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(6); //1 = phi, 2=theta, 3=psi, 4=chi1
		for ( core::Size i(1); i<=3; ++i ) {
			trial_conformations[i][1] =   90;
			trial_conformations[i][2] =   30;
			trial_conformations[i][3] =  -60;
		}
		for ( core::Size i(4); i<=6; ++i ) {
			trial_conformations[i][1] =  -60;
			trial_conformations[i][2] =  -30;
			trial_conformations[i][3] =  -30;
		}
		trial_conformations[ 1][4] = -164.3;
		trial_conformations[ 2][4] =  -56.9;
		trial_conformations[ 3][4] =  128.2;
		trial_conformations[ 4][4] =  -62.5;
		trial_conformations[ 5][4] = -171.0;
		trial_conformations[ 6][4] =   12.6;

		utility::vector1< core::Real > expected_energies( 6 );
		expected_energies[ 1] = -std::log(0.884240);
		expected_energies[ 2] = -std::log(0.115760);
		expected_energies[ 3] = -std::log(1e-6);
		expected_energies[ 4] = -std::log(0.727001);
		expected_energies[ 5] = -std::log(0.269219);
		expected_energies[ 6] = -std::log(0.003781);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tTHETA\tPSI\tCHI1\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_theta( 2, trial_conformations[i][2] );
			pose.set_psi( 2, trial_conformations[i][3] );
			pose.set_chi( 1, 2, trial_conformations[i][4] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}


	/// @brief Confirm that the rotamer wells in the B3N rotamer library really do have
	/// the appropriate score values.
	void test_B3N_rotamer_well_scoring() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[B3N]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 5 > > trial_conformations(20); //1 = phi, 2=theta, 3=psi, 3=chi1, 4=chi2
		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			trial_conformations[i][1] =  -60;
			trial_conformations[i][2] =  -30;
			trial_conformations[i][3] = -150;
		}
		trial_conformations[ 1][4] =  -67.1; trial_conformations[ 1][5] =   -63.2;
		trial_conformations[ 2][4] =  -64.6; trial_conformations[ 2][5] =  -77.2;
		trial_conformations[ 3][4] =  -63.7; trial_conformations[ 3][5] =  112.8;
		trial_conformations[ 4][4] =  -61.9; trial_conformations[ 4][5] =  104.7;
		trial_conformations[ 5][4] =  -61.8; trial_conformations[ 5][5] =  -51.6;
		trial_conformations[ 6][4] =  -71.7; trial_conformations[ 6][5] =  -99.3;
		trial_conformations[ 7][4] =  -67.1; trial_conformations[ 7][5] =   -7.8;
		trial_conformations[ 8][4] =  -76.7; trial_conformations[ 8][5] =    1.2;
		trial_conformations[ 9][4] =  -67.3; trial_conformations[ 9][5] =   80.9;
		trial_conformations[10][4] =  -71.4; trial_conformations[10][5] =   47.8;

		trial_conformations[11][4] =  -65.9; trial_conformations[11][5] =  139.0;
		trial_conformations[12][4] = -169.1; trial_conformations[12][5] = -104.0;
		trial_conformations[13][4] = -173.6; trial_conformations[13][5] =   26.9;
		trial_conformations[14][4] =  178.5; trial_conformations[14][5] =   56.7;
		trial_conformations[15][4] = -173.3; trial_conformations[15][5] =   34.2;
		trial_conformations[16][4] = -154.6; trial_conformations[16][5] =  -75.0;
		trial_conformations[17][4] =  -81.5; trial_conformations[17][5] = -151.4;
		trial_conformations[18][4] =   40.5; trial_conformations[18][5] =  -67.2;
		trial_conformations[19][4] = -158.0; trial_conformations[19][5] =   84.9;
		trial_conformations[20][4] =    7.7; trial_conformations[20][5] =   83.1;

		utility::vector1< core::Real > expected_energies( 20 );
		expected_energies[ 1] = -std::log(0.155934);
		expected_energies[ 2] = -std::log(0.122468);
		expected_energies[ 3] = -std::log(0.098767);
		expected_energies[ 4] = -std::log(0.097278);
		expected_energies[ 5] = -std::log(0.089326);
		expected_energies[ 6] = -std::log(0.079059);
		expected_energies[ 7] = -std::log(0.074263);
		expected_energies[ 8] = -std::log(0.068539);
		expected_energies[ 9] = -std::log(0.060696);
		expected_energies[10] = -std::log(0.041629);

		expected_energies[11] = -std::log(0.030376);
		expected_energies[12] = -std::log(0.028345);
		expected_energies[13] = -std::log(0.019785);
		expected_energies[14] = -std::log(0.017049);
		expected_energies[15] = -std::log(0.005643);
		expected_energies[16] = -std::log(0.004179);
		expected_energies[17] = -std::log(0.003765);
		expected_energies[18] = -std::log(0.001827);
		expected_energies[19] = -std::log(0.000460);
		expected_energies[20] = -std::log(0.000234);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "PHI\tTHETA\tPSI\tCHI1\tCHI2\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_theta( 2, trial_conformations[i][2] );
			pose.set_psi( 2, trial_conformations[i][3] );
			pose.set_chi( 1, 2, trial_conformations[i][4] );
			pose.set_chi( 2, 2, trial_conformations[i][5] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << trial_conformations[i][5] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}


	/// @brief Ensure that rotamer order doesn't matter in a rotamer library.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_equivalence_independent_of_rotamer_order() {
		core::pose::Pose pose_badorder, pose_goodorder;
		core::pose::make_pose_from_sequence(pose_badorder, "GX[DAB]G", "fa_standard");
		core::pose::make_pose_from_sequence(pose_goodorder, "GX[DAB_goodorder]G", "fa_standard");

		// A scorefunction:
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		// Initial, netural backbone conformation:
		for ( core::Size i(1), imax(pose_badorder.total_residue()); i<=imax; ++i ) {
			pose_badorder.set_phi(i, -135);
			pose_badorder.set_phi(i, 135);
			pose_badorder.set_omega(i, 180);
			pose_goodorder.set_phi(i, -135);
			pose_goodorder.set_phi(i, 135);
			pose_goodorder.set_omega(i, 180);
		}

		// The conformations that I'll try:
		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(8); //1 = phi, 2=psi, 3=chi1, 4=chi2
		trial_conformations[1][1] = -61.0; trial_conformations[1][2] = -41.0; trial_conformations[1][3] = 32.0; trial_conformations[1][4] = -53.0;
		trial_conformations[2][1] = -135.0; trial_conformations[2][2] = 135.0; trial_conformations[2][3] = 32.0; trial_conformations[2][4] = -53.0;
		trial_conformations[3][1] = -61.0; trial_conformations[3][2] = -41.0; trial_conformations[3][3] = 72.0; trial_conformations[3][4] = 63.0;
		trial_conformations[4][1] = -61.0; trial_conformations[4][2] = -41.0; trial_conformations[4][3] = -42.0; trial_conformations[4][4] = -122.0;
		trial_conformations[5][1] = -135.0; trial_conformations[5][2] = 135.0; trial_conformations[5][3] = -42.0; trial_conformations[5][4] = -122.0;
		trial_conformations[6][1] = -60.0; trial_conformations[6][2] = -40.0; trial_conformations[6][3] = 59.0; trial_conformations[6][4] = -119.0;
		trial_conformations[7][1] = -60.0; trial_conformations[7][2] = -40.0; trial_conformations[7][3] = 119; trial_conformations[7][4] = -119.0;
		trial_conformations[8][1] = -60.0; trial_conformations[8][2] = -40.0; trial_conformations[8][3] = 30; trial_conformations[8][4] = -80.0;

		utility::vector1< core::Real > energies_badorder, energies_goodorder;

		TR << "\nPHI\tPSI\tCHI1\tCHI2\tBADSCORE\tGOODSCORE\n";
		for ( core::Size i(1); i<=trial_conformations.size(); ++i ) {
			// Set the trial conformation:
			pose_badorder.set_phi(2, trial_conformations[i][1]);
			pose_badorder.set_psi(2, trial_conformations[i][2]);
			pose_badorder.set_chi(1, 2, trial_conformations[i][3]);
			pose_badorder.set_chi(2, 2, trial_conformations[i][4]);
			pose_goodorder.set_phi(2, trial_conformations[i][1]);
			pose_goodorder.set_psi(2, trial_conformations[i][2]);
			pose_goodorder.set_chi(1, 2, trial_conformations[i][3]);
			pose_goodorder.set_chi(2, 2, trial_conformations[i][4]);
			pose_badorder.update_residue_neighbors();
			pose_goodorder.update_residue_neighbors();

			energies_badorder.push_back( sfxn(pose_badorder) );
			energies_goodorder.push_back( sfxn(pose_goodorder) );
			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][4] << "\t" << energies_badorder[i] << "\t" << energies_goodorder[i] << "\n";
		}
		TR << std::endl;

		for ( core::Size i(1); i<=energies_badorder.size(); ++i ) {
			TS_ASSERT_DELTA(energies_badorder[i], energies_goodorder[i], 0.0001);
		}

	}

};
