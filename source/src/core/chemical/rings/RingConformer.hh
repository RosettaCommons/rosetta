// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/RingConformer.hh
/// @brief   Definitions for RingConformer.
/// @author  Labonte <JWLabonte@jhu.edu


#ifndef INCLUDED_core_chemical_rings_RingConformer_HH
#define INCLUDED_core_chemical_rings_RingConformer_HH

// Project header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace chemical {
namespace rings {

/// @brief  A structure for storing information for specific, idealized ring conformers.
struct RingConformer {
	std::string specific_name;  // e.g., "1C4"
	std::string general_name;  // e.g., "chair"

	core::uint degeneracy; // E.g., 1C4 has a degeneracy of 3, since 3CO and 5C2 are equivalent.

	// a list of 1 (for 4-membered rings), 2 (for 5-membered rings), or 3 (for 6-membered rings) Cremer-Pople "ring
	// puckering" parameters
	// TODO: This could be expanded to include parameters for 7-membered rings and larger.
	utility::vector1< core::Real > CP_parameters;  // phi and theta are angles in degrees; q is a distance in Angstroms

	// a list of the 1st n-1 nu angles, where n is the ring size.  Nu angles are internal ring torsions.
	utility::vector1< core::Angle > nu_angles;

	// a list of the n tau angles, where n is the ring size.  Tau angles are internal ring bond angles.
	utility::vector1< core::Angle > tau_angles;
};  // struct RingConformer

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_RingConformer_HH
