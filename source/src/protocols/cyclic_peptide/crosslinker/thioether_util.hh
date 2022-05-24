// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/util.hh
/// @brief Utility functions for setting up thioether cyclization.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_thioether_util_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_thioether_util_hh

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Protocols headers
#include <protocols/simple_moves/DeclareBond.fwd.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

// just from 6qys for now
constexpr core::Real THIOETHER_UTIL_THIOETHER_BOND_LENGTH = 1.827;
constexpr core::Real THIOETHER_UTIL_THIOETHER_BOND_C_ANGLE = 1.781; // 102.028 degrees
constexpr core::Real THIOETHER_UTIL_THIOETHER_BOND_N_ANGLE = 1.960; // 112.284

/// @brief Given a pose and two residues, set up the thioether variant types.
/// @details Sidechainres gets SIDECHAIN_CONJUGATION; ntermres gets ACETYLATED_NTERMINUS_CONNECTION_VARIANT.
/// The pose is modified by this operation.
void
set_up_thioether_variants(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
);

/// @brief Set up the mover that creates thioether lariat bonds.
void
set_up_thioether_bond_mover(
	protocols::simple_moves::DeclareBond & termini,
	core::pose::Pose const & pose,
	core::Size const ntermres,
	core::Size const sidechainres
);

/// @brief Correct the bond angles and bond lenghts for virtual atoms at thioether bonds.
void
correct_thioether_virtuals(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
);

/// @brief Given a pose and two residues to constrain (the first being the residue with the modified N-terminus,
/// and the second being the one with the thiol-containing sidechain), add constraints for a thioether linkage.
void
set_up_thioether_constraints(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
);

} //crosslinker
} //cyclic_peptide
} //protocols


#endif //protocols/cyclic_peptide/crosslinker_thioether_util_hh

