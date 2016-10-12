// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/util.hh
/// @brief   Utility function declarations for scoring carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_util_HH
#define INCLUDED_core_scoring_carbohydrates_util_HH

// Package Headers
#include <core/scoring/carbohydrates/CHIEnergyFunctionLinkageType.hh>
#include <core/scoring/carbohydrates/OmegaPreferenceType.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Header
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace carbohydrates {

/// @brief Get the CHI Energy Function linkage type for phi for a particular residue.
CHIEnergyFunctionLinkageType get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose(
	pose::Pose const & pose,
	core::uint rsd_num );

/// @brief Get the CHI Energy Function linkage type for phi for a particular residue.
CHIEnergyFunctionLinkageType get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose(
	pose::Pose const & pose,
	core::uint rsd_num );

/// @brief Get the omega preference for a particular residue.
OmegaPreferenceType get_omega_preference_for_residue_in_pose(
	pose::Pose const & pose,
	core::uint rsd_num );

/// @brief Get the CHI Energy Function linkage type for the given torsion angle of a particular residue.
CHIEnergyFunctionLinkageType get_CHI_energy_function_linkage_type_for_residue_in_pose(
	id::MainchainTorsionType torsion,
	pose::Pose const & pose,
	core::uint rsd_num );


utility::vector1< CHIEnergyFunctionLinkageType > get_linkage_types_for_dihedral( core::uint torsion_id );

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_carbohydrates_util_HH

