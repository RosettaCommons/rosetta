// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/util.hh
/// @brief Utility functions for setting up lanthionine cyclization.
/// @author Clay Tydings (claiborne.w.tydings@vanderbilt.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_lanthionine_util_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_lanthionine_util_hh

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Protocols headers
#include <protocols/simple_moves/DeclareBond.fwd.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

// from Methionine
constexpr core::Real LANTHIONINE_UTIL_LANTHIONINE_BOND_LENGTH = 1.807;
constexpr core::Real LANTHIONINE_UTIL_LANTHIONINE_BOND_C_ANGLE = 2.011; // 115.230 degrees
constexpr core::Real LANTHIONINE_UTIL_LANTHIONINE_BOND_S_ANGLE = 1.801; // 103.176 degrees
// from MP2/cc-pvtz geometry optimization
constexpr core::Real LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CG_ANGLE = 1.864; // 106.790 degrees SCC to CG
constexpr core::Real LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CA_ANGLE = 1.963; // 112.460 degrees SCC to CA
constexpr core::Real LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CB_ANGLE = 1.927; // 110.440 degrees CA-CB-CG

/// @brief Given a pose and two residues, set up the lanthionine variant types.
/// @details Sidechainres gets SIDECHAIN_CONJUGATION; dalares gets ACETYLATED_NTERMINUS_CONNECTION_VARIANT.
/// The pose is modified by this operation.
void
set_up_lanthionine_variants(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
);

/// @brief Set up the mover that creates lanthionine lariat bonds.
void
set_up_lanthionine_bond_mover(
	protocols::simple_moves::DeclareBond & termini,
	core::pose::Pose const & pose,
	core::Size const dalares,
	core::Size const cysres
);

/// @brief Correct the bond angles and bond lenghts for virtual atoms at lanthionine bonds.
void
correct_lanthionine_virtuals(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
);

/// @brief Given a pose and two residues to constrain (the first being the residue with the modified N-terminus,
/// and the second being the one with the thiol-containing sidechain), add constraints for a lanthionine linkage.
void
set_up_lanthionine_constraints(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
);

} //crosslinker
} //cyclic_peptide
} //protocols


#endif //protocols/cyclic_peptide/crosslinker_lanthionine_util_hh

