// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/min_pack.hh
/// @brief  pack and minimize sidechains
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_min_pack_hh
#define INCLUDED_core_pack_min_pack_hh

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/scmin/AtomTreeCollection.fwd.hh>
#include <core/pack/scmin/SCMinMinimizerMap.fwd.hh>
#include <core/pack/scmin/SidechainStateAssignment.fwd.hh>
#include <core/pack/rotamer_set/ContinuousRotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {

/// @brief Interface function to the minimize-each-rotamer-during-packing sidechain placement algorithm.
void
min_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task,
	bool cartesian=false,
	bool nonideal=false
);

void
min_pack_setup(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	task::PackerTaskOP task,
	bool cartesian,
	bool nonideal,
	rotamer_set::RotamerSetsOP & rotsets,
	scmin::SCMinMinimizerMapOP & scminmap,
	scoring::MinimizationGraphOP & mingraph,
	scmin::AtomTreeCollectionOP & atc,
	optimization::MinimizerOptionsOP & min_options
);

void
min_pack_optimize(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	task::PackerTaskOP task,
	rotamer_set::RotamerSetsOP rotsets,
	scmin::SCMinMinimizerMapOP scminmap,
	scoring::MinimizationGraphOP mingraph,
	scmin::AtomTreeCollectionOP atc,
	optimization::MinimizerOptions const & min_options,
	scmin::SidechainStateAssignment & best_state
);

void
min_pack_place_opt_rotamers_on_pose(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	rotamer_set::RotamerSetsOP rotsets,
	scmin::AtomTreeCollectionOP atc,
	scmin::SidechainStateAssignment const & best_state,
	Real start_score
);

/// @brief Interface to a version of the packer that uses very little memory
/// and simultaneously is able to go off rotamer and explore more of sidechain
/// conformation space.  Quite a bit faster than the min-packer.
void
off_rotamer_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task
);

void
off_rotamer_pack_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task,
	rotamer_set::ContinuousRotamerSetsOP & rotsets,
	scmin::AtomTreeCollectionOP & atc,
	interaction_graph::SimpleInteractionGraphOP & ig
);

void
off_rotamer_pack_optimize(
	rotamer_set::ContinuousRotamerSets const & rotsets,
	scmin::AtomTreeCollectionOP atc,
	interaction_graph::SimpleInteractionGraph & ig,
	scmin::SidechainStateAssignment & best_state
);

void
off_rotamer_pack_update_pose(
	pose::Pose & pose,
	rotamer_set::ContinuousRotamerSets const & rotsets,
	scmin::AtomTreeCollectionOP atc,
	scmin::SidechainStateAssignment const & best_state
);

} // namespace pack
} // namespace core

#endif
