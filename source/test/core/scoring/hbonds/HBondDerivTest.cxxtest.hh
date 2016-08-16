// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBondDerivTest.cxxtest.hh
/// @brief  Test evaluation of hbond potentials for each HBEvalType across parameter ranges.
/// @author Matthew O'Meara


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Package headers
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>

//Project headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


//Numeric headers
#include <numeric/constants.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>

// Utility headers
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>


using namespace core;
using namespace conformation;
using namespace scoring;
using namespace hbonds;

using pose::Pose;
using conformation::Residue;


static basic::Tracer TR("core.scoring.hbonds.HBondDerivTest.cxxtest");

class HBondDerivTest : public CxxTest::TestSuite {

public:
	void setUp() {
		core_init();

	}


	void tearDown(){}

	void test_dummy_test() {}

	void fail_test_chi_0_pole() {
		do_test_chi_0_pole(hbacc_AHX);
		do_test_chi_0_pole(hbacc_CXL);
	}

	void do_test_chi_0_pole(HBAccChemType acc_chem_type) {
		HBondDatabaseCOP database(HBondDatabase::get_database());
		HBondOptions hbond_options;
		hbond_options.measure_sp3acc_BAH_from_hvy(true);

		bool apply_chi_torsion_penalty;
		HBGeoDimType AHD_geometric_dimension;
		// Get Analytic derivative
		Real energy, dE_dr, dE_dxD, dE_dxH;

		Real min_energy(100), max_energy(-100);

		HBEvalTuple hbt( hbdon_GDE, acc_chem_type, seq_sep_other );
		Real AHdis(2.0);
		Real xD(0.99999);
		Real xH(0.2);

		Size n_steps = 100;
		for ( Size i=0; i < n_steps; ++i ) {
			Real chi(numeric::constants::d::pi_2 * i / n_steps);

			hbond_compute_energy(
				*database,
				hbond_options,
				hbt, AHdis, xD, xH, chi,
				energy,
				apply_chi_torsion_penalty,
				AHD_geometric_dimension,
				dE_dr, dE_dxD, dE_dxH);

			min_energy = std::min(min_energy, energy);
			max_energy = std::max(max_energy, energy);
		}

		TR
			<< "For Acceptor Chemical Type: " << HBondTypeManager::name_from_acc_chem_type(acc_chem_type) << std::endl
			<< "Circling the BAH=0 pole at angle of acos(1-.999) leads to an energy difference of " << max_energy - min_energy << std::endl;
		TS_ASSERT(max_energy - min_energy < .0001);
	}

	void do_not_test_hbond_deriv_a(){
		HBondDatabaseCOP database(HBondDatabase::get_database());
		HBondOptions hbond_options;


		// A few bonds that have been causing some trouble...
		{
			HBEvalTuple hbt( hbdon_GDE, hbacc_AHX, seq_sep_other );
			Real AHdis(2.0626650237011339);
			Real xD(0.81221037050947442);
			Real xH(0.37311327017559648);
			Real chi(0.0);
			do_hbond_deriv_test( database, hbond_options, hbt, AHdis, xD, xH, chi );
		}

		{
			HBEvalTuple hbt( hbdon_GDE, hbacc_HXL, seq_sep_other );
			Real AHdis(2.722940506144047);
			Real xD(0.62920743488168385);
			Real xH(0.15382334407518322);
			Real chi(0.0);
			do_hbond_deriv_test( database, hbond_options,  hbt, AHdis, xD, xH, chi );
		}

		{
			HBEvalTuple hbt( hbdon_PBA, hbacc_PBA, seq_sep_other );
			Real AHdis(1.7);
			Real xD(0.98);
			Real xH(-0.48);
			Real chi(0.0);
			do_hbond_deriv_test( database, hbond_options,  hbt, AHdis, xD, xH, chi );
		}

	}

	void test_hbond_sp2_deriv_no_fade_energy()
	{
		using namespace core::scoring::hbonds;
		//To debug at full precision set:
		//std::cout.precision(16);
		HBondDatabaseCOP database(HBondDatabase::get_database("sp2_params"));


		HBondOptions hboptions;
		hboptions.use_sp2_chi_penalty(true);

		HBEvalTuple hbt = HBEvalTuple( hbdon_PBA, hbacc_PBA, seq_sep_other );
		Real xD(1.0-.02);
		for ( Real AHdis = 1.7; AHdis <= 1.9; AHdis += .1 ) {
			for ( Real xH = MAX_xH - 0.02; xH > MIN_xH + 0.02; xH -= .3 ) {
				for ( Real chi = 0.0; chi < 2*numeric::constants::d::pi; chi += 0.3 ) {
					do_hbond_deriv_test(database, hboptions, hbt, AHdis, xD, xH, chi);
				}
			}
		}

	}

	void test_hbond_sp3acc_deriv_no_fade_energy()
	{
		using namespace core::scoring::hbonds;
		//To debug at full precision set:
		//std::cout.precision(16);
		HBondDatabaseCOP database(HBondDatabase::get_database("sp2_params"));


		HBondOptions hboptions;
		hboptions.use_sp2_chi_penalty(true);

		HBEvalTuple hbt = HBEvalTuple( hbdon_HXL, hbacc_HXL, seq_sep_other );
		Real xD(1.0-.02);
		for ( Real AHdis = 1.7; AHdis <= 1.9; AHdis += .1 ) {
			for ( Real xH = MAX_xH - 0.02; xH > MIN_xH + 0.02; xH -= .3 ) {
				for ( Real chi = 0.0; chi < 2*numeric::constants::d::pi; chi += 0.3 ) {
					do_hbond_deriv_test(database, hboptions, hbt, AHdis, xD, xH, chi);
				}
			}
		}

	}

	void do_no_test_hbond_sp2_deriv_fade_energy()
	{
		using namespace core::scoring::hbonds;
		//To debug at full precision set:
		//std::cout.precision(16);
		HBondDatabaseCOP database(HBondDatabase::get_database("OLF_params_11a"));


		HBondOptions hboptions;
		hboptions.use_sp2_chi_penalty(true);
		hboptions.fade_energy(true);

		HBEvalTuple hbt = HBEvalTuple( hbdon_PBA, hbacc_PBA, seq_sep_other );
		Real xD(1.0-.02);
		for ( Real AHdis = 1.7; AHdis <= 1.9; AHdis += .1 ) {
			for ( Real xH = MAX_xH - 0.02; xH > MIN_xH + 0.02; xH -= .3 ) {
				for ( Real chi = 0.0; chi < 2*numeric::constants::d::pi; chi += 0.5 ) {
					do_hbond_deriv_test(database, hboptions, hbt, AHdis, xD, xH, chi, .00005, 5, false, .001);
				}
			}
		}

	}


	// This is a numeric_deriv type test: Compare the derivative at a
	// point to the slope of planes constructed to approximate the
	// tangent at the point.

	void do_hbond_deriv_test(
		HBondDatabaseCOP database,
		HBondOptions hbond_options,
		HBEvalTuple hbt,
		Real const AHdis, // acceptor proton distance
		Real const xD, // -cos(180-theta), where theta is defined by Tanja K.
		Real const xH, // cos(180-phi), where phi is defined by Tanja K.
		Real const chi, // AB2-AB-A-H dihdral angle for sp2 hybridized acceptors
		Real const increment = 0.0005,
		Size const n_increment = 5,
		bool verbose = false,
		Real deriv_tolerance = .001){

		bool analytic_apply_chi_torsion_penalty;
		bool numeric_apply_chi_torsion_penalty;
		HBGeoDimType AHD_geometric_dimension;
		// Get Analytic derivative
		Real energy, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi;
		hbond_compute_energy(
			*database,
			hbond_options,
			hbt, AHdis, xD, xH, chi,
			energy,
			analytic_apply_chi_torsion_penalty,
			AHD_geometric_dimension,
			dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);

		// Only check derivatives that have a hope of becoming an hbond
		if ( energy > .1 ) return;

		Real BAH(numeric::constants::d::pi - acos(xH));
		Real e_low, e_high, dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH, dummy_dE_dBAH, dummy_dE_dchi;
		Real numeric_deriv, deriv, deriv_dev;

		if ( verbose ) {
			if ( !analytic_apply_chi_torsion_penalty ) {
				TR
					<< "Coordinates:"
					<< "AHdis: " << AHdis << ", "
					<< "xD: " << xD << ", "
					<< "xH: " << xH << ", "
					<< "BAH: " << BAH << ", "
					<< "chi: " << chi << " (!apply_chi_torsion_penalty)" << std::endl;
			} else {
				TR
					<< "Coordinates:"
					<< "AHdis: " << AHdis << ", "
					<< "xD: " << xD << ", "
					<< "xH: " << xH << ", "
					<< "BAH: " << BAH << ", "
					<< "chi: " << chi << " (apply_chi_torsion_penalty)" << std::endl;
			}
			TR << "Energy: " << energy << std::endl;
		}

		// Check AHdis
		Real test_AHdis;
		deriv_dev = 10000.0;
		deriv = dE_dr;
		for ( Size j = 1, factor=1; j <= n_increment; ++j ) {
			factor*=2;
			test_AHdis = AHdis - factor * increment;
			hbond_compute_energy(
				*database,
				hbond_options,
				hbt, test_AHdis, xD, xH, chi,
				e_low,
				numeric_apply_chi_torsion_penalty,
				AHD_geometric_dimension,
				dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);

			test_AHdis = AHdis + factor * increment;
			hbond_compute_energy(
				*database,
				hbond_options,
				hbt, test_AHdis, xD, xH, chi,
				e_high,
				numeric_apply_chi_torsion_penalty,
				AHD_geometric_dimension,
				dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);
			numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
			deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
			Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
				deriv / numeric_deriv );
			if ( verbose &&
					( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
				TR
					<< "dim:AHdis"
					<< "\tstep:" << factor*increment
					<< "\tanalyitic_deriv: " << deriv
					<< "\tnumeric_deriv:" << numeric_deriv
					<< "\te_low:" << e_low
					<< "\te_high:" << e_high
					<< "\tratio:" << ratio << std::endl;
			}
		}

		if ( deriv_dev > deriv_tolerance ) {
			TR << "DERIV ERROR: dim:AHdis derivative deviation:" << deriv_dev << std::endl;
			TS_ASSERT(false);
		}

		if ( analytic_apply_chi_torsion_penalty ) {
			// Check BAH
			Real test_xH;
			deriv_dev = 10000.0;
			deriv = dE_dxH + dE_dBAH / sqrt( 1 - xH * xH);
			for ( Size j = 1, factor=1; j <= n_increment; ++j ) {

				factor*=2;
				test_xH = xH - factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, test_xH, chi,
					e_low,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH, dummy_dE_dBAH, dummy_dE_dchi);

				test_xH = xH + factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, test_xH, chi,
					e_high,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH, dummy_dE_dBAH, dummy_dE_dchi);
				numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
				deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
				Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
					deriv / numeric_deriv );
				if ( verbose &&
						( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
					TR << "dim:BAH"
						<< "\tstep:" << factor*increment
						<< "\tanalyitic_deriv: " << deriv
						<< "\tnumeric_deriv:" << numeric_deriv
						<< "\te_low:" << e_low
						<< "\te_high:" << e_high
						<< "\tratio:" << ratio << std::endl;
				}
			}

			if ( deriv_dev > deriv_tolerance ) {
				TR << "DERIV ERROR: dim:BAH derivative deviation:" << deriv_dev << std::endl;
				TS_ASSERT(false);
			}

			// Check chi
			Real test_chi;
			deriv_dev = 10000.0;
			deriv = dE_dchi;
			for ( Size j = 1, factor=1; j <= n_increment; ++j ) {

				factor*=2;
				test_chi = chi - factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, xH, test_chi,
					e_low,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH, dummy_dE_dBAH, dummy_dE_dchi);

				test_chi = chi + factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, xH, test_chi,
					e_high,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH, dummy_dE_dBAH, dummy_dE_dchi);
				numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
				deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
				Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
					deriv / numeric_deriv );
				if ( verbose &&
						( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
					TR << "dim:chi"
						<< "\tstep:" << factor*increment
						<< "\tanalyitic_deriv: " << deriv
						<< "\tnumeric_deriv:" << numeric_deriv
						<< "\te_low:" << e_low
						<< "\te_high:" << e_high
						<< "\tratio:" << ratio << std::endl;
				}
			}

			if ( deriv_dev > deriv_tolerance ) {
				TR << "DERIV ERROR: dim:chi derivative deviation:" << deriv_dev << std::endl;
				TS_ASSERT(false);
			}


		} else {

			// Check cosBAH
			Real test_xH;
			deriv_dev = 10000.0;
			deriv = dE_dxH;
			for ( Size j = 1, factor=1; j <= n_increment; ++j ) {
				factor*=2;
				test_xH = xH - factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, test_xH, chi,
					e_low,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);

				test_xH = xH + factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, xD, test_xH, chi,
					e_high,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);
				numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
				deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
				Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
					deriv / numeric_deriv );
				if ( verbose &&
						( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
					TR
						<< "dim:cosBAH"
						<< "\tstep:" << factor*increment
						<< "\tanalyitic_deriv: " << deriv
						<< "\tnumeric_deriv:" << numeric_deriv
						<< "\te_low:" << e_low
						<< "\te_high:" << e_high
						<< "\tratio:" << ratio << std::endl;
				}
			}

			if ( deriv_dev > deriv_tolerance ) {
				TR << "DERIV ERROR: dim:cosBAH derivative deviation:" << deriv_dev << std::endl;
				TS_ASSERT(false);
			}
		}

		if ( AHD_geometric_dimension == hbgd_cosAHD ) {
			// Check cosAHD
			Real test_xD;
			deriv_dev = 10000.0;
			deriv = dE_dxD;
			for ( Size j = 1, factor=1; j <= n_increment; ++j ) {
				factor*=2;
				test_xD = xD - factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, test_xD, xH, chi,
					e_low,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);

				test_xD = xD + factor * increment;
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, test_xD, xH, chi,
					e_high,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);
				numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
				deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
				Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
					deriv / numeric_deriv );
				if ( verbose &&
						( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
					TR
						<< "dim:cosAHD"
						<< "\tstep:" << factor*increment
						<< "\tanalyitic_deriv: " << deriv
						<< "\tnumeric_deriv:" << numeric_deriv
						<< "\te_low:" << e_low
						<< "\te_high:" << e_high
						<< "\tratio:" << ratio << std::endl;
				}
			}

			if ( deriv_dev > deriv_tolerance ) {
				TR << "DERIV ERROR: dim:cosAHD derivative deviation:" << deriv_dev << std::endl;
				TS_ASSERT(false);
			}
		} else if ( AHD_geometric_dimension == hbgd_AHD ) {
			// Check AHD
			Real test_xD;
			Real test_AHD;
			deriv_dev = 10000.0;
			Real AHD( numeric::constants::d::pi - acos(xD) );
			deriv = dE_dxD; // the dE_dxD value is actually dE/dAHD
			for ( Size j = 1, factor=1; j <= n_increment; ++j ) {
				factor*=2;
				test_AHD = AHD - factor * increment;
				test_xD = cos(numeric::constants::d::pi - test_AHD);
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, test_xD, xH, chi,
					e_low,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);

				test_AHD = AHD + factor * increment;
				test_xD = cos(numeric::constants::d::pi - test_AHD);
				hbond_compute_energy(
					*database,
					hbond_options,
					hbt, AHdis, test_xD, xH, chi,
					e_high,
					numeric_apply_chi_torsion_penalty,
					AHD_geometric_dimension,
					dummy_dE_dr, dummy_dE_dxD, dummy_dE_dxH);
				numeric_deriv = ( e_high - e_low ) / ( factor * 2 * increment );
				deriv_dev = std::min( deriv_dev, std::abs( deriv - numeric_deriv ) );
				Real const ratio( std::abs( numeric_deriv ) < .001 ? 0.0 :
					deriv / numeric_deriv );
				if ( verbose &&
						( std::abs(numeric_deriv) > 0.001 || std::abs(deriv) > 0.001 ) ) {
					TR
						<< "dim:AHD"
						<< "\tstep:" << factor*increment
						<< "\tanalyitic_deriv: " << deriv
						<< "\tnumeric_deriv:" << numeric_deriv
						<< "\te_low:" << e_low
						<< "\te_high:" << e_high
						<< "\tratio:" << ratio << std::endl;
				}
			}

			if ( deriv_dev > deriv_tolerance ) {
				TR << "DERIV ERROR: dim:cosAHD derivative deviation:" << deriv_dev << std::endl;
				TS_ASSERT(false);
			}
		}
	}

	void do_not_test_f1f2_deriv(){
		HBondDatabaseCOP database(HBondDatabase::get_database());
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunctionOP score_function( get_score_function() );
		score_function->score( pose );

		HBondSet hbond_set( pose.n_residue() );

		fill_hbond_set(pose,
			true,
			hbond_set);

		for ( Size i=1; i <= hbond_set.nhbonds(); ++i ) {
			TR << "HBond " << hbond_set.hbond(i) << std::endl;
			do_f1f2_deriv_test(database, hbond_set.hbond_options(), pose, hbond_set.hbond(i));
			TR << std::endl;
		}

	}


	// In the Abe and Go derivative calculation, each score term must
	// return two vectors, f1 and f2.  To understand what the f1 and f2
	// vectors are, imagine the hydrogen and donor atoms are fixed and
	// space and the acceptor and base atoms are rotated an infantesimal
	// dphi about an arbitrary vector Eab:=(Vb-Va).  How does the energy
	// of the hydrogen bond change?  This is dE/dphi the quantity we are
	// interested in.
	//

	// dE/dphi = Eab . f1  +  Eab x Vb . f2

	void
	do_f1f2_deriv_test(
		HBondDatabaseCOP database,
		HBondOptions const & hbond_options,
		Pose const & pose,
		HBond const & hbond
	)
	{
		using namespace numeric::deriv;
		using std::abs;
		using std::cos;
		using std::sqrt;
		using numeric::constants::d::pi;

		// Compute reference f1 f2 derivative vectors;
		Vector f1( hbond.derivs().h_deriv.f1() );
		Vector f2( hbond.derivs().h_deriv.f2() );


		// Compute the three dE_d<hbond_geo_dim>
		// Copied from hb_energy_deriv
		const Residue acc_res( pose.residue(hbond.acc_res()));
		const Vector Axyz(acc_res.atom(hbond.acc_atm()).xyz());
		const Vector Bxyz(acc_res.atom(acc_res.atom_base(hbond.acc_atm())).xyz());
		const Vector B2xyz(acc_res.atom(acc_res.abase2(hbond.acc_atm())).xyz());
		const Residue don_res( pose.residue(hbond.don_res()));
		const Vector Hxyz(don_res.atom(hbond.don_hatm()).xyz());
		const Vector HDunit(create_don_orientation_vector(don_res, hbond.don_hatm()));
		chemical::Hybridization acc_hybrid( get_hbe_acc_hybrid( hbond.eval_type() ) );
		Vector BAunit;
		Vector PBxyz;
		make_hbBasetoAcc_unitvector( hbond_options, acc_hybrid, Axyz, Bxyz, B2xyz, PBxyz, BAunit );


		// compute hbond geometric dimensions
		Vector AH;
		AH = Hxyz - Axyz;

		Real const AHdis2( AH.length_squared() );

		Real const AHdis = sqrt(AHdis2);
		Real const inv_AHdis = 1.0f / AHdis;
		Vector AHunit(AH * inv_AHdis);

		Real const xD =            // cos(180-theta) = cos(thetaD)
			dot( AHunit, HDunit );

		Real const xH =            // cos(180-psi) = cos(thetaH)
			dot( BAunit, AHunit );

		Real chi( 0 );
		if ( hbond_options.use_sp2_chi_penalty() &&
				get_hbe_acc_hybrid( hbond.eval_type() ) == chemical::SP2_HYBRID &&
				B2xyz != Vector(-1.0, -1.0, -1.0) ) {
			chi = numeric::dihedral_radians( Hxyz, Axyz, Bxyz, B2xyz );
		} else if ( hbond.eval_tuple().acc_type() == hbacc_AHX || hbond.eval_tuple().acc_type() == hbacc_HXL ) {
			/// Bxyz really is the heavy atom base and B2xyz really is the hydroxyl hydrogen
			/// this is guaranteed by the hbond_measure_sp3acc_BAH_from_hvy flag.
			chi = numeric::dihedral_radians( Hxyz, Axyz, Bxyz, B2xyz );
		}


		Real dE_dxH, dE_dxD, dE_dr, dummy_energy;
		bool apply_chi_torsion_penalty;
		HBGeoDimType AHD_geometric_dimension;
		hbond_compute_energy(
			*database,
			hbond_options,
			hbond.eval_tuple(),
			AHdis, xD, xH, chi,
			dummy_energy,
			apply_chi_torsion_penalty,
			AHD_geometric_dimension,
			dE_dr, dE_dxD, dE_dxH);


		// Compute f1 f2 using generic computations
		// It isn't done this way in the HBond code because it is slower

		// The temp_AHdis, temp_xD and temp_xH are only within 10^-3 of
		// AHdis, xD and xH because the numeric::deriv versions adjust the
		// value when angles get close to being degenerate

		Vector new_f1(0.0), new_f2(0.0);
		Vector temp_f1, temp_f2;

		Real temp_AHdis;
		distance_f1_f2_deriv(Hxyz, Axyz, temp_AHdis, temp_f1, temp_f2);
		TS_ASSERT( abs( temp_AHdis - AHdis ) < .001 );
		new_f1 += temp_f1*dE_dr;
		new_f2 += temp_f2*dE_dr;

		Vector new_AHdis_f1( temp_f1 *dE_dr );
		Vector new_AHdis_f2( temp_f2 *dE_dr );

		TR << "new_AHdis_f1: " << new_AHdis_f1[0] << " "<< new_AHdis_f1[1] << " " << new_AHdis_f1[2] << " new_AHdis_f2:" << new_AHdis_f2[0] << " " << new_AHdis_f2[1] << " " << new_AHdis_f2[2] << std::endl;

		Real theta;
		angle_p1_deriv(  Axyz, Hxyz, Hxyz-HDunit, theta, temp_f1, temp_f2);

		if ( AHD_geometric_dimension == hbgd_cosAHD ) {

			Real comp_xD(- cos(pi - theta ) );
			TS_ASSERT( abs( xD - comp_xD ) < 0.001);

			new_f1 += temp_f1*dE_dxD*sin(theta);
			new_f2 += temp_f2*dE_dxD*sin(theta);

			Vector new_xD_f1( temp_f1*dE_dxD*sin(theta));
			Vector new_xD_f2( temp_f2*dE_dxD*sin(theta));

			TR << "new_xD_f1: " << new_xD_f1[0] << " "<< new_xD_f1[1] << " " << new_xD_f1[2] << " new_xD_f2:" << new_xD_f2[0] << " " << new_xD_f2[1] << " " << new_xD_f2[2] << std::endl;
		} else if ( AHD_geometric_dimension == hbgd_AHD ) {

			TS_ASSERT( abs( xD - theta ) < 0.001);

			new_f1 += temp_f1*dE_dxD;
			new_f2 += temp_f2*dE_dxD;

			Vector new_xD_f1( temp_f1*dE_dxD);
			Vector new_xD_f2( temp_f2*dE_dxD);

			TR << "new_xD_f1: " << new_xD_f1[0] << " "<< new_xD_f1[1] << " " << new_xD_f1[2] << " new_xD_f2:" << new_xD_f2[0] << " " << new_xD_f2[1] << " " << new_xD_f2[2] << std::endl;
			TR << "||f1-new_xD_f1||: " << (f1-new_xD_f1).length() << std::endl;
			TR << "||f2-new_xD_f2||: " << (f2-new_xD_f2).length() << std::endl;
			TS_ASSERT( (f1-new_xD_f1).length() < .000000001 );
			TS_ASSERT( (f2-new_xD_f2).length() < .000000001 );


		} else {
			utility_exit_with_message("Invalid AHD_geometric_dimension");
		}

		Real phi;
		angle_p1_deriv( Hxyz, Axyz,
			Axyz - BAunit, //virtual atom base
			phi,
			temp_f1, temp_f2);
		Real comp_xH( std::cos( pi - phi ) );
		TS_ASSERT( std::abs(xH - comp_xH) < .001 );

		new_f1 += temp_f1*dE_dxH*sin(phi);
		new_f2 += temp_f2*dE_dxH*sin(phi);


		Vector new_xH_f1( temp_f1*dE_dxH*sin(phi));
		Vector new_xH_f2( temp_f2*dE_dxH*sin(phi));

		TR << "new_xH_f1: " << new_xH_f1[0] << " "<< new_xH_f1[1] << " " << new_xH_f1[2] << " new_xH_f2:" << new_xH_f2[0] << " " << new_xH_f2[1] << " " << new_xH_f2[2] << std::endl;


		TR << "f1: " << f1[0] << " "<< f1[1] << " " << f1[2] << " f2:" << f2[0] << " " << f2[1] << " " << f2[2] << std::endl;
		TR << "new_f1: " << new_f1[0] << " " << new_f1[1] << " " << new_f1[2] << " new_f2:" << new_f2[0] << " " << new_f2[1] << " " << new_f2[2] << std::endl;
		TR << "||f1-new_f1||: " << (f1-new_f1).length() << std::endl;
		TR << "||f2-new_f2||: " << (f2-new_f2).length() << std::endl;
		TS_ASSERT( (f1-new_f1).length() < .000000001 );
		TS_ASSERT( (f2-new_f2).length() < .000000001 );
	}

};
