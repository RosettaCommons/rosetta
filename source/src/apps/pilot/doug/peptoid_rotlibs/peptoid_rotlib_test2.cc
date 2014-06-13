// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/doug/peptoid_rotlibs/peptoid_rotlib_test2.cc
/// @brief Test the peptoid rotlibs
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// core headers
#include <core/init.hh>
#include <core/types.hh>

// numeric functions
#include <numeric/numeric.functions.hh>
#include <numeric/angle.functions.hh>
#include <basic/basic.hh>

// c++ headers
#include <iostream>
#include <string>

// namespaces
using namespace core;

const Real OMG_BINRANGE(10.00);
const Real PHI_BINRANGE(10.00);
const Real PSI_BINRANGE(10.00);
const Size N_OMG_BINS(7);
const Size N_PHI_BINS(36);
const Size N_PSI_BINS(36);
const Real CIS_OMG_LOWER_RANGE(-30.0);
const Real CIS_OMG_UPPER_RANGE(30.0);
const Real TRANS_OMG_LOWER_RANGE(150.0);
const Real TRANS_OMG_UPPER_RANGE(-150.0);
const Real CIS_UPPER_TRANS_LOWER_SPLIT(90.0);
const Real TRANS_UPPER_CIS_LOWER_SPLIT(-90.0);

void
bin_angle(
  Real const angle_start,
  Real const angle_step,
  Real const angle_range,
  Size const nbins,
  Real const ang,
  Size & bin_lower,
  Size & bin_upper,
  Real & angle_alpha
);

void
get_omgphipsi_bins(
	Real omg,
	Real phi,
	Real psi,
	Size & omgbin,
	Size & phibin,
	Size & psibin,
	Size & omgbin_next,
	Size & phibin_next,
	Size & psibin_next,
	Real & omg_alpha,
	Real & phi_alpha,
	Real & psi_alpha
);

int
main( int argc, char * argv [] )
{
	// init options, rng, etc.
	core::init(argc, argv);

	Size omg_bin_lower, omg_bin_upper, phi_bin_lower, phi_bin_upper, psi_bin_lower, psi_bin_upper;
	Real omg, phi, psi, omg_alpha, phi_alpha, psi_alpha;

	omg_bin_lower = omg_bin_upper = phi_bin_lower = phi_bin_upper = psi_bin_lower = psi_bin_upper = 0;
	omg = phi = psi = omg_alpha = phi_alpha = psi_alpha = 0;

	for ( omg = -180.00; omg <= 180; omg += 1 ) {
	//for ( phi = -180.00; phi <= 180; phi += 5 ) {
		get_omgphipsi_bins( omg, phi, psi,
			omg_bin_lower, phi_bin_lower, psi_bin_lower,
			omg_bin_upper, phi_bin_upper, psi_bin_upper,
			omg_alpha, phi_alpha, psi_alpha );

		//std::cout << omg << "\t" << omg_bin_lower << "\t" << omg_bin_upper << "\t" << omg_alpha << "\t"
		//					<< phi << "\t" << phi_bin_lower << "\t" << phi_bin_upper << "\t" << phi_alpha << "\t"
		//						<< psi << "\t" << psi_bin_lower << "\t" << psi_bin_upper << "\t" << psi_alpha
		//					<< std::endl;
	}

	return 0;
}

void
bin_angle(
  Real const angle_start,
  Real const angle_step,
  Real const angle_range,
  Size const nbins,
  Real const ang,
  Size & bin_lower,
  Size & bin_upper,
  Real & angle_alpha
)
{
  /// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
  /// though it is supposed to return values in the range [-180, 180).
	assert( angle_start <= ang && ang <= angle_start + angle_range );
	assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

  Real real_bin_lower = ( ang - angle_start ) / angle_step;
  Size bin_prev = static_cast< Size > ( real_bin_lower );
  bin_lower = 1 + numeric::mod( bin_prev, nbins );
  bin_upper = numeric::mod( bin_lower, nbins ) + 1;
  angle_alpha = ( (ang - angle_start ) - ( bin_prev * angle_step ) ) / angle_step;
}

void
get_omgphipsi_bins(
	Real omg,
	Real phi,
	Real psi,
	Size & omgbin,
	Size & phibin,
	Size & psibin,
	Size & omgbin_next,
	Size & phibin_next,
	Size & psibin_next,
	Real & omg_alpha,
	Real & phi_alpha,
	Real & psi_alpha
)
{

	/*

		Need something that works for the following ranges...

		[   30,   90 ) use cis_range_upper
		[   90,  150 ) use trans_range_lower
		[  150,  180 ) data
		[ -180, -150 ) data
		[ -150,  -90 ) use trans_range_upper
		[  -90,  -30 ) use cis_range_lower
		[  -30,    0 ) data
		[    0,   30 ) data


	 */

	Real nnpad_omg( numeric::nonnegative_principal_angle_degrees(omg) );

	if ( omg >= 0.00 && omg < CIS_OMG_UPPER_RANGE ) {
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 1 << std::endl;
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= CIS_OMG_UPPER_RANGE && omg < CIS_UPPER_TRANS_LOWER_SPLIT ) { // [   30,   90 ) use cis_range_upper
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 1 << std::endl;
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(CIS_OMG_UPPER_RANGE), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= CIS_UPPER_TRANS_LOWER_SPLIT && omg < TRANS_OMG_LOWER_RANGE ) { // [   90,  150 ) use trans_range_lower
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 2 << std::endl;
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(TRANS_OMG_LOWER_RANGE), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= TRANS_OMG_LOWER_RANGE && omg < 180.00 ) { //	[  150, 180 ) use omg
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 3 << std::endl;
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );

	} else if (omg >= -180.00 && omg < TRANS_OMG_UPPER_RANGE ) { // [  180, -150 use nnpad
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 3.5 << std::endl;
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= TRANS_OMG_UPPER_RANGE && omg < TRANS_UPPER_CIS_LOWER_SPLIT ) { // [ -150,  -90 ) use trans_range_upper
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 4 << std::endl;
		bin_angle( 150.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::nonnegative_principal_angle_degrees(TRANS_OMG_UPPER_RANGE), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= TRANS_UPPER_CIS_LOWER_SPLIT && omg < CIS_OMG_LOWER_RANGE) { // [  -90,  -30 ) use cis_range_lower
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 5 << std::endl;
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(CIS_OMG_LOWER_RANGE), omgbin, omgbin_next, omg_alpha );

	} else if ( omg >= CIS_OMG_LOWER_RANGE && omg < CIS_OMG_UPPER_RANGE ) { // [  -30,    30 ) use omg
		//std::cout << omg << "\t" << nnpad_omg << "\t" << 6 << std::endl;
		bin_angle( -30.0, OMG_BINRANGE, 70.0, N_OMG_BINS, numeric::principal_angle_degrees(omg), omgbin, omgbin_next, omg_alpha );
	}


	bin_angle( -180.0, PHI_BINRANGE, 360.0, N_PHI_BINS, basic::periodic_range( phi, 360 ), phibin, phibin_next, phi_alpha );
	bin_angle( -180.0, PSI_BINRANGE, 360.0, N_PSI_BINS, basic::periodic_range( psi, 360 ), psibin, psibin_next, psi_alpha );

	std::cout << omg << "\t" << omgbin << "\t" << omgbin_next << "\t" << omg_alpha << "\t"
						<< phi << "\t" << phibin << "\t" << phibin_next << "\t" << phi_alpha << "\t"
						<< psi << "\t" << psibin << "\t" << psibin_next << "\t" << psi_alpha
						<< std::endl;
}
