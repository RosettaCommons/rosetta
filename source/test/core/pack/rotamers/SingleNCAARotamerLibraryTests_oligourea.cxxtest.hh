// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/rotamers/SingleNCAARotamerLibraryTests_oligourea.cxxtest.hh
/// @brief  Unit tests for the non-canonical rotamer libraries, with oligoureas.
/// @details Split into a separate unit test suite because of runtime.
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

static basic::Tracer TR("SingleNCAARotamerLibraryTests_oligourea");


class SingleNCAARotamerLibraryTests_oligourea : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		// This is just leucine, but treated as though it's a non-canonical amino acid.
		core_init_with_additional_options("-extra_res_fa core/pack/rotamers/LEU_NCAA.params core/pack/rotamers/DAB_corrected_order.params");
	}

	void tearDown(){

	}

	/// @brief Confirm that the rotamer wells in the oligourea-valine rotamer library really do have
	/// the appropriate score values.
	void test_OU3_VAL_rotamer_well_scoring() {
		core::pose::Pose pose, dpose;
		core::pose::make_pose_from_sequence(pose, "GX[OU3_VAL]G", "fa_standard");
		core::pose::make_pose_from_sequence(dpose, "GX[DOU3_VAL]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 4 > > trial_conformations(6); //1 = phi, 2=theta, 3=psi, 4=chi1
		for ( core::Size i(1); i<=3; ++i ) {
			trial_conformations[i][1] =  -60;
			trial_conformations[i][2] =  -50;
			trial_conformations[i][3] =   30;
		}
		for ( core::Size i(4); i<=6; ++i ) {
			trial_conformations[i][1] =  110;
			trial_conformations[i][2] =  160;
			trial_conformations[i][3] =  -50;
		}

		trial_conformations[ 1][4] = 171.9;
		trial_conformations[ 2][4] = -62.8;
		trial_conformations[ 3][4] =  61.4;
		trial_conformations[ 4][4] = 175.2;
		trial_conformations[ 5][4] =  85.2;
		trial_conformations[ 6][4] = -62.8;

		utility::vector1< core::Real > expected_energies( 6 );
		expected_energies[ 1] = -std::log(0.999391);
		expected_energies[ 2] = -std::log(0.000480);
		expected_energies[ 3] = -std::log(0.000129);
		expected_energies[ 4] = -std::log(0.735716);
		expected_energies[ 5] = -std::log(0.264284);
		expected_energies[ 6] = -std::log(1e-6);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "Testing L-version:" << std::endl;
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

		TR << "Testing D-version:" << std::endl;
		TR << "PHI\tTHETA\tPSI\tCHI1\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			dpose.set_phi( 2, -1.0*trial_conformations[i][1] );
			dpose.set_theta( 2, -1.0*trial_conformations[i][2] );
			dpose.set_psi( 2, -1.0*trial_conformations[i][3] );
			dpose.set_chi( 1, 2, -1.0*trial_conformations[i][4] );

			core::Real const actual_energy( sfxn(dpose) );

			TR << -1.0*trial_conformations[i][1] << "\t" << -1.0*trial_conformations[i][2] << "\t" << -1.0*trial_conformations[i][3] << "\t" << -1.0*trial_conformations[i][4] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

	/// @brief Confirm that the rotamer wells in the oligourea-proline rotamer library really do have
	/// the appropriate score values.
	void test_OU3_PRO_rotamer_well_scoring() {
		core::pose::Pose pose, dpose;
		core::pose::make_pose_from_sequence(pose, "GX[OU3_PRO]G", "fa_standard");
		core::pose::make_pose_from_sequence(dpose, "GX[DOU3_PRO]G", "fa_standard");

		utility::vector1< utility::fixedsizearray1< core::Real, 6 > > trial_conformations(6); //1 = phi, 2=theta, 3=psi, 4=chi1, 5=chi2, 6=chi3
		for ( core::Size i(1); i<=2; ++i ) {
			trial_conformations[i][1] =  -60;
			trial_conformations[i][2] = -110;
			trial_conformations[i][3] = -120;
		}
		for ( core::Size i(3); i<=4; ++i ) {
			trial_conformations[i][1] =  -50;
			trial_conformations[i][2] = -170;
			trial_conformations[i][3] = -180;
		}
		for ( core::Size i(5); i<=6; ++i ) {
			trial_conformations[i][1] = -140;
			trial_conformations[i][2] =    0;
			trial_conformations[i][3] =  -40;
		}

		trial_conformations[ 1][4] =  28.8;
		trial_conformations[ 2][4] = -25.1;
		trial_conformations[ 3][4] = -23.9;
		trial_conformations[ 4][4] =  27.1;
		trial_conformations[ 5][4] = -25.1;
		trial_conformations[ 6][4] =  27.6;

		trial_conformations[ 1][5] = -35.8;
		trial_conformations[ 2][5] =  36.4;
		trial_conformations[ 3][5] =  34.3;
		trial_conformations[ 4][5] = -34.8;
		trial_conformations[ 5][5] =  36.3;
		trial_conformations[ 6][5] = -34.5;

		trial_conformations[ 1][6] =  27.6;
		trial_conformations[ 2][6] = -31.9;
		trial_conformations[ 3][6] = -30.3;
		trial_conformations[ 4][6] =  27.9;
		trial_conformations[ 5][6] = -31.8;
		trial_conformations[ 6][6] =  27.0;

		utility::vector1< core::Real > expected_energies( 6 );
		expected_energies[ 1] = -std::log(0.738847);
		expected_energies[ 2] = -std::log(0.261153);
		expected_energies[ 3] = -std::log(0.552691);
		expected_energies[ 4] = -std::log(0.447309);
		expected_energies[ 5] = -std::log(0.832961);
		expected_energies[ 6] = -std::log(0.167039);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_dun, 1.0 );

		TR << "Testing L-OU3_PRO:" << std::endl;
		TR << "PHI\tTHETA\tPSI\tCHI1\tCHI2\tCHI3\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			pose.set_phi( 2, trial_conformations[i][1] );
			pose.set_theta( 2, trial_conformations[i][2] );
			pose.set_psi( 2, trial_conformations[i][3] );
			pose.set_chi( 1, 2, trial_conformations[i][4] );
			pose.set_chi( 2, 2, trial_conformations[i][5] );
			pose.set_chi( 3, 2, trial_conformations[i][6] );

			core::Real const actual_energy( sfxn(pose) );

			TR << trial_conformations[i][1] << "\t" << trial_conformations[i][2] << "\t" << trial_conformations[i][3] << "\t" << trial_conformations[i][5] << "\t" << trial_conformations[i][6] << "\t" << trial_conformations[i][4] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;

		TR << "Testing D-OU3_PRO:" << std::endl;
		TR << "PHI\tTHETA\tPSI\tCHI1\tCHI2\tCHI3\tEXPECT\tACTUAL" << std::endl;

		for ( core::Size i(1), imax(trial_conformations.size()); i<=imax; ++i ) {
			dpose.set_phi( 2, -trial_conformations[i][1] );
			dpose.set_theta( 2, -trial_conformations[i][2] );
			dpose.set_psi( 2, -trial_conformations[i][3] );
			dpose.set_chi( 1, 2, -trial_conformations[i][4] );
			dpose.set_chi( 2, 2, -trial_conformations[i][5] );
			dpose.set_chi( 3, 2, -trial_conformations[i][6] );

			core::Real const actual_energy( sfxn(dpose) );

			TR << -trial_conformations[i][1] << "\t" << -trial_conformations[i][2] << "\t" << -trial_conformations[i][3] << "\t" << -trial_conformations[i][5] << "\t" << -trial_conformations[i][6] << "\t" << -trial_conformations[i][4] << "\t" << expected_energies[i] << "\t" << actual_energy << std::endl;
			TS_ASSERT_DELTA( expected_energies[i], actual_energy, 1e-6 );
		}

		TR << std::endl;
	}

};
