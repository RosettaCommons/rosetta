// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/RamachandranEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::RamachandranEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>

//Auto Headers
#include <utility/vector1.hh>

static basic::Tracer TR( "DunbrackEnergyTests" );

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;

class DunbrackEnergyTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void do_not_test_dunbrack_deriv_check()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( fa_dun,   0.4 );
		kinematics::MoveMap movemap; movemap.set_bb( true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void bin_angle_mod(
		Real const angle_start,
		Real const angle_step,
		Real const ASSERT_ONLY( angle_range ),
		Size const nbins,
		Real const ang,
		Size & bin_lower,
		Size & bin_upper,
		Real & angle_alpha
	) const {
		/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
		/// though it is supposed to return values in the range [-180, 180).
		debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
		// temp -- by merge-time we won't call this function at all.
		//debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

		Real real_bin_lower = ( ang - angle_start ) / angle_step;
		TR << "rea_bin_lower " << real_bin_lower << std::endl;
		Size bin_prev = static_cast< Size > ( real_bin_lower );
		TR << "bin_prev      " << bin_prev << std::endl;
		bin_lower = 1 + numeric::mod( bin_prev, nbins );
		TR << "bin_lower     " << bin_lower << std::endl;
		bin_upper = numeric::mod( bin_lower, nbins ) + 1;
		TR << "bin_upper     " << bin_upper << std::endl;
		angle_alpha = ( (ang - angle_start ) - ( bin_prev * angle_step ) ) / angle_step;
	}

	void
	get_bb_bins_mod(
		Real & bbs,
		Size & bb_bin,
		Size & bb_bin_next,
		Real & bb_alpha
	) const
	{

		bbs = numeric::principal_angle_degrees( bbs );

		if ( bbs >= 0.00 && bbs < 30 ) {
			bin_angle_mod( -30.0, 10, 70.0, 14, numeric::principal_angle_degrees( bbs ), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 3; bb_bin_next += 3;
		} else if ( bbs >= 30 && bbs < 90 ) { // [   30,   90 ) use cis_range_upper
			bin_angle_mod( -30.0, 10, 70.0, 14, numeric::principal_angle_degrees(30), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 3; bb_bin_next += 3;
		} else if ( bbs >= 90 && bbs < 150 ) { // [   90,  150 ) use trans_range_lower
			bin_angle_mod( 150.0, 10, 70.0, 14, numeric::principal_angle_degrees(150), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 10; bb_bin_next += 10;
		} else if ( bbs >= 150 && bbs <= 180.00 ) { // [  150, 180 ] use bbs
			bin_angle_mod( 150.0, 10, 70.0, 14, numeric::nonnegative_principal_angle_degrees( bbs ), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 10; bb_bin_next += 10;
		} else if ( bbs >= -180.00 && bbs < -150 ) { // [  180, -150 ) use nnpad
			bin_angle_mod( 150.0, 10, 70.0, 14, numeric::nonnegative_principal_angle_degrees( bbs ), bb_bin, bb_bin_next, bb_alpha );
			bb_bin -= 4; bb_bin_next -= 4;
		} else if ( bbs >= -150 && bbs < -90 ) { // [ -150,  -90 ) use trans_range_upper
			bin_angle_mod( 150.0, 10, 70.0, 14, numeric::nonnegative_principal_angle_degrees(-150), bb_bin, bb_bin_next, bb_alpha );
			bb_bin -= 4; bb_bin_next -= 4;
		} else if ( bbs >= -90 && bbs < -30 ) { // [  -90,  -30 ) use cis_range_lower
			bin_angle_mod( -30.0, 10, 70.0, 14, numeric::principal_angle_degrees(-30), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 3; bb_bin_next += 3;
		} else if ( bbs >= -30 && bbs < 30 ) { // [  -30,    30 ) use bbs
			bin_angle_mod( -30.0, 10, 70.0, 14, numeric::principal_angle_degrees( bbs ), bb_bin, bb_bin_next, bb_alpha );
			bb_bin += 3; bb_bin_next += 3;
		}
		if ( bb_bin == 0 ) bb_bin = 14;
		if ( bb_bin > 14 ) bb_bin -= 14;
		if ( bb_bin_next > 14 ) bb_bin_next -= 14;
	}

	void test_get_bb_bins_range()
	{

		Size bb_bin( 999999 );
		Size bb_bin_next( 999999 );
		Real bin_alpha( 999999 );

		for ( Real i(-180); i <= 180; i += 2 ) {
			get_bb_bins_mod( i, bb_bin, bb_bin_next, bin_alpha );
			TR << "BLARG angle: "  << std::setw(4) << i << " bb_bin " << bb_bin << " bb_bin_next " << bb_bin_next << " bin_alpha " << bin_alpha << std::endl;
		}
		TR << "DOUG DOUG DOUG STOP" << std::endl;
	}

	void test_peptoid_energy_snap()
	{
		// A peptoid should basically 'nearest-neighbor interpolate' when
		// caught between cis and trans.
		Pose peptoid_cis;
		make_pose_from_sequence(peptoid_cis, "A[ALA:NtermProteinFull]X[601]A[ALA:CtermProteinFull]", core::chemical::FA_STANDARD, true);

		ScoreFunction sfxn;
		sfxn.set_weight( fa_dun, 0.75 );

		for ( Size ii = 1; ii <= peptoid_cis.size(); ++ii ) {
			peptoid_cis.set_phi( ii, -150 );
			peptoid_cis.set_psi( ii, -150 );
			peptoid_cis.set_omega( ii, 180 );
			peptoid_cis.set_chi( 1, 2, 180 );
			peptoid_cis.set_chi( 2, 2, 180 );
		}
		Pose peptoid_cislike = peptoid_cis;
		peptoid_cis.set_omega( 1, 30 );
		peptoid_cislike.set_omega( 1, 40 );

		sfxn( peptoid_cis );
		sfxn( peptoid_cislike );

		//TR << " 30 vs 40 " << std::endl;
		TS_ASSERT_DELTA( peptoid_cis.energies().total_energy(), peptoid_cislike.energies().total_energy(), 1e-6 );

		peptoid_cis.set_omega( 1, 150 );
		peptoid_cislike.set_omega( 1, 140 );

		sfxn( peptoid_cis );
		sfxn( peptoid_cislike );

		//TR << "150 vs 140 " << std::endl;
		TS_ASSERT_DELTA( peptoid_cis.energies().total_energy(), peptoid_cislike.energies().total_energy(), 1e-6 );

		peptoid_cis.set_omega( 1, -30 );
		peptoid_cislike.set_omega( 1, -40 );

		sfxn( peptoid_cis );
		sfxn( peptoid_cislike );

		//TR << "-30 vs -40 " << std::endl;
		TS_ASSERT_DELTA( peptoid_cis.energies().total_energy(), peptoid_cislike.energies().total_energy(), 1e-6 );

		peptoid_cis.set_omega( 1, -150 );
		peptoid_cislike.set_omega( 1, -140 );

		sfxn( peptoid_cis );
		sfxn( peptoid_cislike );

		//TR << "-150 vs -140 " << std::endl;
		TS_ASSERT_DELTA( peptoid_cis.energies().total_energy(), peptoid_cislike.energies().total_energy(), 1e-6 );
	}

};
