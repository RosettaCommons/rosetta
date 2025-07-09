// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/util.hh
/// @brief  Utilities for ligand docking
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_ligand_docking_util_hh
#define INCLUDED_protocols_ligand_docking_util_hh

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

#include <string>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief Move the center of specified chain to the desired_centroid
void move_ligand_to_desired_centroid(
	std::string const & chain,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
);

/// @brief Move the center of specified multiple chains to the desired_centroid
void move_ligand_to_desired_centroid(
	utility::vector1<std::string> const & chains,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
);

/// @brief Move the center of the object(s) downstream of jump_id to the desired_centroid
void move_ligand_to_desired_centroid(
	core::Size const jump_id,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
);

/// @brief Move the neighbor atom of the specified multiple chains to the desired_position
void move_ligand_neighbor_to_desired_position(
	utility::vector1<std::string> const & chains,
	core::Vector const & desired_position,
	core::pose::Pose & pose
);

/// @brief Move the neighbor atom of the specified chain to the desired_position
void move_ligand_neighbor_to_desired_position(
	std::string const & chain,
	core::Vector const & desired_position,
	core::pose::Pose & pose
);

/// @brief Annotate a string/real pair.
/// Will place it in the JD2 job, and if JD2 is not active, will also place it in the pose extra scores.
void
add_string_real_pair( core::pose::Pose & pose, std::string const & string, core::Real real );

} //namespace ligand_docking
} //namespace protocols

#endif
