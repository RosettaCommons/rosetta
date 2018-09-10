// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/util.hh
/// @brief  Utility functions for ApproximateBuriedUnsatPenalty
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_approximate_buried_unsat_penalty_util_hh
#define INCLUDED_core_pack_guidance_scoreterms_approximate_buried_unsat_penalty_util_hh

// Unit headers

// Package headers
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/pack_basic/RotamerSetsBase.hh>


// #include <basic/datacache/CacheableStringFloatMathMatrixMap.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {



scoring::hbonds::graph::AtomLevelHBondGraphOP
hbond_graph_from_partial_rotsets(
	pose::Pose const & pose_in,
	pack::rotamer_set::RotamerSetsOP const & original_rotsets,
	scoring::ScoreFunctionOP const & scorefxn,
	pack::rotamer_set::RotamerSetsOP & complete_rotsets_out,
	utility::vector1<bool> & position_had_rotset,
	float minimum_hb_cut = 0
);

}
}
}
}

#endif
