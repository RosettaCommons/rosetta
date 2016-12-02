// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/OmegaPreferenceType.hh
/// @brief   Enumerator definition for OmegaPreferenceType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_OmegaPreferenceType_HH
#define INCLUDED_core_scoring_carbohydrates_OmegaPreferenceType_HH

namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief    Labels for the omega torsion preference type used to select the proper OmegaPreferencesFunction form.
/// @details  In saccharide residues where the hydroxyl group of the carbon atom two carbons previous to the exocyclic
/// carbon is equatorial, the "gauche effect" occurs, in which the exocyclic torsion angle prefers one of the two gauche
/// orientations instead of the expected anti configuration.\n
enum OmegaPreferenceType {
	ANTI = 1,       // The gauche effect is not in place, i.e., the anti rotamer is preferred.
	GAUCHE_EFFECT,  // The gauche effect is in place, i.e., the gauche rotamers are preferred.
	
	N_OMEGA_PREFERENCE_TYPES = GAUCHE_EFFECT,
	
	PREFERENCE_NA   // The gauche effect is not applicable here.
};

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_carbohydrates_OmegaPreferenceType_HH
