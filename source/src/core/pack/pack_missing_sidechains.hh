// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/pack_missing_sidechains.hh
/// @brief  header for subroutine to run rotamer trials on residues missing sidechain density in a PDB
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_core_pack_pack_missing_sidechains_hh
#define INCLUDED_core_pack_pack_missing_sidechains_hh

// Project headers
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {

/// @brief This function runs rotamer trials on residues missing sidechain density (as described by the AtomID_Mask)
void
pack_missing_sidechains(
	pose::Pose & pose,
	id::AtomID_Mask const& missing
);

/// @brief return vector of bools with true for each residue that has >=1 atom in to_repack that is not VIRTUAL, ORBS or LPbb
bool figure_out_repackable_residues(
  core::pose::Pose & pose,
	core::id::AtomID_Mask const& to_repack,
	utility::vector1_bool& repackable
);


} //namespace pack
} //namespace core

#endif //INCLUDED_core_pack_pack_missing_sidechains_HH
