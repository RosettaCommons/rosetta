// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/InteractionGraphFactory.hh
/// @brief  Interation graph factory class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_InteractionGraphFactory_hh
#define INCLUDED_core_pack_interaction_graph_InteractionGraphFactory_hh

// Unit headers
#include <core/pack/interaction_graph/InteractionGraphFactory.fwd.hh>

// Package headers
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class InteractionGraphFactory {
public:

	/// @brief Create appropiate InteractionGraph instance for given
	/// packer task, rotamer set, and pose.
	static
	InteractionGraphBaseOP
	create_interaction_graph(
		task::PackerTask const & packer_task,
		rotamer_set::RotamerSets const & rotsets,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::Graph const & packer_neighbor_graph
	);

	/// @brief Create and initialize two-body interaction graph for the given
	/// pose, rotamer sets, packer task and score function.
	///
	/// Call only valid after:
	/// Pose has been scored by scorefxn.
	/// Pose residue neighbors updated.
	/// ScoreFunction setup for packing.
	/// Rotamer sets built.
	static
	InteractionGraphBaseOP
	create_and_initialize_two_body_interaction_graph(
		task::PackerTask const & packer_task,
		rotamer_set::RotamerSets & rotsets,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph);

	/// @brief Create and initialize annealable graph for the given
	/// pose, rotamer sets, packer task and score function. Initalizes
	/// two-body interaction graph, as well as other annealable graphs
	/// sepecified by the given score function and task.
	/// @details
	/// Call only valid after:
	/// Pose has been scored by scorefxn.
	/// Pose residue neighbors updated.
	/// ScoreFunction setup for packing.
	/// Rotamer sets built.
	/// @note Pose is nonconst as there may still be data to cache in
	/// the pose at this point.
	static
	AnnealableGraphBaseOP
	create_and_initialize_annealing_graph(
		task::PackerTask const & packer_task,
		rotamer_set::RotamerSets & rotsets,
		pose::Pose & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph);

private:

	/// @brief Clear specific types of cached information in the pose that the ResidueArrayAnnealableEneriges use.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	static
	void
	clear_cached_residuearrayannealableenergy_information(
		core::pose::Pose &pose
	);

	/// @brief Apply IGEdgeReweights to IG if specified in packer task.
	static
	void
	setup_IG_res_res_weights(
		pose::Pose const & pose,
		task::PackerTask const & task,
		rotamer_set::RotamerSets const & rotsets,
		interaction_graph::InteractionGraphBase & ig
	);
};

}
}
}

#endif
