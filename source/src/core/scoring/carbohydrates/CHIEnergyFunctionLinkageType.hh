// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/CHIEnergyFunctionLinkageType.hh
/// @brief   Enumerator definition for CHIEnergyFunctionLinkageType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_CHIEnergyFunctionLinkageType_HH
#define INCLUDED_core_scoring_carbohydrates_CHIEnergyFunctionLinkageType_HH

namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief    Labels for the CHI Energy Function linkage type of the carbohydrate phi or psi angle.
/// @details  The CHI energy functions depend on the type of linkage.\n
/// Parameters for the Gaussian functions that compose the CHI Energy Functions are stored in a vector indexed by these
/// labels.
enum CHIEnergyFunctionLinkageType {
	FIRST_LINK_TYPE = 1,
	ALPHA_LINKS = 1,     // used to describe alpha linkages (for scoring phi)
	BETA_LINKS,          // used to describe beta linkages (for scoring phi)
	_2AX_3EQ_4AX_LINKS,  // used to describe ->2-axial, ->3-equatorial, or ->4-axial linkages (for scoring psi)
	_2EQ_3AX_4EQ_LINKS,  // used to describe ->2-equatorial, ->3-axial, or ->4-equatorial linkages (for scoring psi)
	ALPHA6_LINKS,        // used to describe alpha exocyclic linkages (for scoring psi)
	BETA6_LINKS,         // used to describe beta exocyclic linkages (for scoring psi)
	N_LINK_TYPES = BETA6_LINKS,

	LINKAGE_NA           // a null for searching
};

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_carbohydrates_CHIEnergyFunctionLinkageType_HH
