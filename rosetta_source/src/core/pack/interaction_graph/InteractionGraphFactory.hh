// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/InteractionGraphFactory.hh
/// @brief  Interation graph factory class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_InteractionGraphFactory_hh
#define INCLUDED_core_pack_interaction_graph_InteractionGraphFactory_hh

// Unit headers
#include <core/pack/interaction_graph/InteractionGraphFactory.fwd.hh>

// Package headers
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class InteractionGraphFactory {
public:

	static
	InteractionGraphBaseOP
	create_interaction_graph(
		task::PackerTask const &,
		rotamer_set::RotamerSets const &,
		pose::Pose const &,
		scoring::ScoreFunction const &
	);
};

}
}
}

#endif
