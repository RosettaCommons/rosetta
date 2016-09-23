// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/SwitchResidueTypeSet.hh
/// @brief Functions for switching the residue type set of a pose
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_util_SwitchResidueTypeSet_hh
#define INCLUDED_core_util_SwitchResidueTypeSet_hh

// Unit headers

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace util {

/// @details the function allows a pose to use a different residue_type_set to
/// represent all its residues, such as from fullatom residues to centroid
/// residues, or vice versa. During the switch, corresponding atoms will be
/// copied. Redundant atoms will be removed (in case from fullatom to centroid)
/// and missing atoms will be built by ideal geometry (in the case from centroid
/// to fullatom).
void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & type_set,
	bool allow_sloppy_match = false
);

/// @details the function allows a pose to use a different residue_type_set to
/// represent all its residues, such as from fullatom residues to centroid
/// residues, or vice versa. During the switch, corresponding atoms will be
/// copied. Redundant atoms will be removed (in case from fullatom to centroid)
/// and missing atoms will be built by ideal geometry (in the case from centroid
/// to fullatom).
void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	core::chemical::TypeSetCategory type_set_type,
	bool allow_sloppy_match = false
);

/// @details the function allows a pose to use a different residue_type_set to
/// represent all its residues, such as from fullatom residues to centroid
/// residues, or vice versa. During the switch, corresponding atoms will be
/// copied. Redundant atoms will be removed (in case from fullatom to centroid)
/// and missing atoms will be built by ideal geometry (in the case from centroid
/// to fullatom).
void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	std::string const & type_set_name,
	bool allow_sloppy_match = false
);

//////////////
// Functions primarily for internal use
// (Mainly to cut down on the giant if-then-else function)
////////////

void
switch_to_centroid_rot_set(
	core::pose::Pose & pose,
	core::conformation::symmetry::SymmetryInfoCOP symm_info,
	core::chemical::ResidueTypeSet const & rsd_set,
	bool allow_sloppy_match = false
);

/// @brief Rebuild disulfides after a transition to a full atom ResidueTypeSet
void
rebuild_fa_disulfides(
	core::pose::Pose & pose
);

} // util
} // core

#endif //INCLUDED_core_util_switchresiduetypeset_HH

