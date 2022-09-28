// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/util.hh
/// @brief Utilities for DrugDesign
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_util_hh
#define INCLUDED_protocols_drug_design_util_hh

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>

namespace protocols {
namespace drug_design {

/// @brief Reassign the coordinates of `residue` to best match that of `target` using the provided mapping
void align_residues(
	core::conformation::Residue & residue,
	core::conformation::Residue const & target,
	core::chemical::IndexIndexMapping const & map // from target to residue
);

/// @brief Replace the residue at the given position in the pose with the new residue type, using the atom correspondence in the mapping
/// The method used can be specified ("default") being whatever is currently thought of as "best".
void
place_new_restype(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype,
	core::chemical::IndexIndexMapping const & map, // from original residuetype in pose to new_restype
	std::string const & method = "default"
);

/// @brief Replaces the residue at the given position in the pose with the new residue type.
/// Does no alignment - assumes that the internal coordinates of the new restype are properly set.
void
place_new_restype_no_align(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype
);

/// @brief Replaces the residue at the given position in the pose with the new residue type,
/// using the rotamer with the best alignment of mapped atoms and a fixed rigid body position of the neighbor atom
void
place_new_restype_rotamer_align(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype,
	core::chemical::IndexIndexMapping const & map // from original residuetype in pose to new_restype
);


} //protocols
} //drug_design


#endif //protocols/drug_design_util_hh

