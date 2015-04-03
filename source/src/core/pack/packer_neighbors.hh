// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/packer_neighbors.hh
/// @brief  creates a graph that describes the possible connectivity induced by designing-in larger side chains
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_packer_neighbors_hh
#define INCLUDED_core_pack_packer_neighbors_hh


// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/types.hh>

// Utility Headers

#include <utility/vector1.hh>


namespace core {
namespace pack {

graph::GraphOP
create_packer_graph(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task
);

graph::GraphOP
create_packer_graph(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	core::Size total_nodes,
	utility::vector1< Distance > const & residue_radii
);

utility::vector1< Distance >
find_residue_max_radii(
	pose::Pose const & pose,
	task::PackerTaskCOP task
);

void
pack_scorefxn_pose_handshake(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn
);


} // namespace core
} // namespace pack

#endif
