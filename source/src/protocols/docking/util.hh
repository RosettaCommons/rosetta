// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util
/// @brief protocols that are specific to docking low resolution
/// @details
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_docking_util_hh
#define INCLUDED_protocols_docking_util_hh

// Package headers
#include <protocols/docking/types.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>


namespace protocols {
namespace docking {

/// @brief Setup FoldTree for docking across an interface.
/// The partners are described by a string for the partner chains (using pdb_chain identification) separated by "_".
/// The foldtree is set up such that the jump points are at the center of masses of the two partners.
void
setup_foldtree(
	core::pose::Pose & pose,
	std::string const & partner_chainID,
	DockJumps & movable_jumps );

/// @brief Setup foldtree for docking across an interface.
/// The partners are described by a vector of booleans of the same length as the number of residues in the pose.
/// The FoldTree is set up such that the jump points are at the center of masses of the two partners.
void
setup_foldtree(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & partner1,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft);

} // docking
} // protocols

#endif
