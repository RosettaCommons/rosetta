// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/GreekDistance.hh
/// @brief   Enumeration definitions for GreekDistance.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_GreekDistance_HH
#define INCLUDED_core_chemical_GreekDistance_HH

namespace core {
namespace chemical {

/// @brief  Enumerators for the Greek distance from the atom with the functional group of highest priority.
enum GreekDistance {
	PRIMARY_ATOM = 0,
	ALPHA_ATOM,
	BETA_ATOM,
	GAMMA_ATOM,
	DELTA_ATOM,
	EPSILON_ATOM,
	ZETA_ATOM,
	ETA_ATOM,
	THETA_ATOM,
	IOTA_ATOM,
	KAPPA_ATOM,
	LAMBDA_ATOM,
	MU_ATOM,
	NU_ATOM,
	XI_ATOM,
	OMICRON_ATOM,
	PI_ATOM,
	SIGMA_ATOM,
	TAU_ATOM,
	UPSILON_ATOM,
	PHI_ATOM,
	CHI_ATOM,
	PSI_ATOM,
	NA_GREEK_DISTANCE = 1023
};

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_GreekDistance_HH
