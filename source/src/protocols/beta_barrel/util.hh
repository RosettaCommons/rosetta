// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/util.hh
/// @brief  Utility functions for beta-barrel construction.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_util_hh
#define INCLUDED_protocols_beta_barrel_util_hh

#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace beta_barrel {

/// @brief Compute per-strand azimuthal position, axial stagger, and inversion flag
/// for each strand in a beta-barrel, given global barrel parameters.
/// @param[in]  n_strands       Number of strands in the barrel.
/// @param[in]  shear_number    Shear number of the barrel.
/// @param[in]  z1              Rise per residue along the strand axis (Angstroms).
/// @param[in]  antiparallel    If true, even-numbered strands are inverted.
/// @param[out] delta_omega0_values  Azimuthal position of each strand (radians).
/// @param[out] delta_z0_values      Axial offset of each strand (Angstroms).
/// @param[out] invert_values        Inversion flag for each strand.
void compute_strand_positions(
	core::Size n_strands,
	core::Size shear_number,
	core::Real z1,
	bool antiparallel,
	utility::vector1< core::Real > & delta_omega0_values,
	utility::vector1< core::Real > & delta_z0_values,
	utility::vector1< bool > & invert_values
);

} //namespace beta_barrel
} //namespace protocols

#endif //INCLUDED_protocols_beta_barrel_util_hh
