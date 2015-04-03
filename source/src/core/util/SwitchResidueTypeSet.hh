// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/util/SwitchResidueTypeSet.hh
/// @brief Functions for switching the residue type set of a pose
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_util_SwitchResidueTypeSet_hh
#define INCLUDED_core_util_SwitchResidueTypeSet_hh

// Unit headers

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
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
	std::string const & type_set_name,
	bool allow_sloppy_match = false
	);

} // util
} // core

#endif //INCLUDED_core_util_switchresiduetypeset_HH

