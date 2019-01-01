// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/parameters.hh
/// @brief   utility functions that return constant values of the dipolar coupling constant
///          or RDC scaling factor for the choice of the RDC atom type
/// @details last Modified: 11/24/18
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Utility headers
#include <core/types.hh>
#include <core/scoring/nmr/types.hh>

#ifndef INCLUDED_core_scoring_nmr_rdc_parameters_HH
#define INCLUDED_core_scoring_nmr_rdc_parameters_HH

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

/// @brief returns the ideal bond length, in case this is not calculated from the pose
inline
Real
rdc_ideal_bond_length(RDC_TYPE const & type) {
	Real blen(0);
	switch ( type ) {
	case RDC_TYPE_NH :
		blen = 1.041;
		break;
	case RDC_TYPE_NCO :
		blen = 1.329;
		break;
	case RDC_TYPE_NCA :
		blen = 1.458;
		break;
	case RDC_TYPE_CAHA :
		blen = 1.107;
		break;
	case RDC_TYPE_CAHN :
		blen = 2.560;
		break;
	case RDC_TYPE_COHN :
		blen = 2.107;
		break;
	case RDC_TYPE_CACO :
		blen = 1.525;
		break;
	case RDC_TYPE_CACB :
		blen = 1.525;
		break;
	}
	return blen;
}

/// @brief returns the ratio (gA * gB) / (r * r * r) relative to N-H
///        this can be used to scale all RDC values relative to NH values
///        in which case one general form of the RDC equation can be used
///        Note that the 15N gyromagnetic ratio is treated as positive here.
///        To account for the different signs of the NH and CH dipolar coupling
///        constant, the option correct_sign can be used which changes the sign of Dconst.
inline
Real
rdc_scaling_factor_toNH(RDC_TYPE const & type) {
	Real sfac(1.0);
	switch ( type ) {
	case RDC_TYPE_NH :
		sfac = 1.00000;
		break;
	case RDC_TYPE_NCO :
		sfac = 0.12084;
		break;
	case RDC_TYPE_NCA :
		sfac = 0.09152;
		break;
	case RDC_TYPE_CAHA :
		sfac = 2.06278;
		break;
	case RDC_TYPE_CAHN :
		sfac = 0.16679;
		break;
	case RDC_TYPE_COHN :
		sfac = 0.29916;
		break;
	case RDC_TYPE_CACO :
		sfac = 0.19839;
		break;
	case RDC_TYPE_CACB :
		sfac = 0.19839;
		break;
	}
	return sfac;
}

/// @brief returns the ratio (gA * gB) / (r * r * r) relative to CA-HA
///        this can be used to scale all RDC values relative to CAHA values
///        in which case one general form of the RDC equation can be used
///        Note that the 15N gyromagnetic ratio is treated as positive here.
///        To account for the different signs of the NH and CH dipolar coupling
///        constant, the option correct_sign can be used which changes the sign of Dconst.
inline
Real
rdc_scaling_factor_toCH(RDC_TYPE const & type) {
	Real sfac(1.0);
	switch ( type ) {
	case RDC_TYPE_NH :
		sfac = 0.48478;
		break;
	case RDC_TYPE_NCO :
		sfac = 0.05858;
		break;
	case RDC_TYPE_NCA :
		sfac = 0.04437;
		break;
	case RDC_TYPE_CAHA :
		sfac = 1.00000;
		break;
	case RDC_TYPE_CAHN :
		sfac = 0.08086;
		break;
	case RDC_TYPE_COHN :
		sfac = 0.14503;
		break;
	case RDC_TYPE_CACO :
		sfac = 0.09617;
		break;
	case RDC_TYPE_CACB :
		sfac = 0.09617;
		break;
	}
	return sfac;
}

/// @brief returns the dipolar coupling constant for two spins A and B in Hz/Ang^3
///        D_const = -(gA * gB * mu0 * hbar)/(4 * pi * pi)
///        where gA, gB are the gyromagnetic ratios of the two spins and mu0 and  hbar
///        are the vacuum permeability and Planck's constant / 2*pi
///        From: Kramer et al. (2004) Concepts Magn Reson Part A, 21, 10-21
///        It is a common NMR convention to treat the sign of the 15N gyromagnetic ratio
///        as positive. As a consequence the sign of dipolar couplings of experiments
///        not involving 15N is flipped. To account for that, I flip here the sign of the
///        CH and CC dipolar coupling constants. To turn this behavior off, the user can set
///        the option "correct_sign" to true.
inline
Real
rdc_D_const(
	RDC_TYPE const & type,
	bool correct_sign = false
) {
	Real D_AB(1.0);
	switch ( type ) {
	case RDC_TYPE_NH :
		D_AB = 24349.878;
		break;
	case RDC_TYPE_NCO :
		D_AB = 6122.400;
		break;
	case RDC_TYPE_NCA :
		D_AB = 6122.400;
		break;
	case RDC_TYPE_CAHA :
		D_AB = correct_sign ? -60400.558 : 60400.558;
		break;
	case RDC_TYPE_CAHN :
		D_AB = correct_sign ? -60400.558 : 60400.558;
		break;
	case RDC_TYPE_COHN :
		D_AB = correct_sign ? -60400.558 : 60400.558;
		break;
	case RDC_TYPE_CACO :
		D_AB = correct_sign ? -15186.785 : 15186.785;
		break;
	case RDC_TYPE_CACB :
		D_AB = correct_sign ? -15186.785 : 15186.785;
		break;
	}
	return D_AB;
}


/// @brief returns the maximal dipolar coupling constant (in Hz) for a given spin pair A and B
///        taking into account a predefined distance r between A and B
///        D_max = -(gA * gB * mu0 * hbar)/(4 * pi * pi * r^3)
///        where gA, gB are the gyromagnetic ratios of the two spins and mu0 and  hbar
///        are the vacuum permeability and Planck's constant / 2*pi
///        From: Kramer et al. (2004) Concepts Magn Reson Part A, 21, 10-21
///        It is a common NMR convention to treat the sign of the 15N gyromagnetic ratio
///        as positive. As a consequence the sign of dipolar couplings of experiments
///        not involving 15N is flipped. To account for that, I flip here the sign of the
///        CH and CC dipolar coupling constants. To turn this behavior off, the user can set
///        the option "correct_sign" to true.
inline
Real
rdc_D_max(
	RDC_TYPE const & type,
	bool correct_sign = false
) {
	Real D_AB(1.0);
	switch ( type ) {
	case RDC_TYPE_NH : // r = 1.041 Ang.
		D_AB = 21584.630;
		break;
	case RDC_TYPE_NCO : // r = 1.329 Ang.
		D_AB = 2608.235;
		break;
	case RDC_TYPE_NCA : // r = 1.458 Ang.
		D_AB = 1975.373;
		break;
	case RDC_TYPE_CAHA : // r = 1.107 Ang.
		D_AB = correct_sign ? -44524.401 : 44524.401;
		break;
	case RDC_TYPE_CAHN : // r = 2.560 Ang.
		D_AB = correct_sign ? -3600.154 : 3600.154;
		break;
	case RDC_TYPE_COHN : // r = 2.107 Ang.
		D_AB = correct_sign ? -6457.246 : 6457.246;
		break;
	case RDC_TYPE_CACO : // r = 1.525 Ang.
		D_AB = correct_sign ? -4282.095 : 4282.095;
		break;
	case RDC_TYPE_CACB : // r = 1.525 Ang.
		D_AB = correct_sign ? -4282.095 : 4282.095;
		break;
	}
	return D_AB;
}

} // rdc
} // nmr
} // scoring
} // core

#endif /* INCLUDED_core_scoring_nmr_rdc_parameters_HH */
