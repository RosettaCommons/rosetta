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
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack_basic/RotamerSetsBase.hh>

#include <basic/datacache/CacheableUint64MathMatrixFloatMap.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {


basic::datacache::CacheableUint64MathMatrixFloatMapOP
three_body_approximate_buried_unsat_calculation(
	pose::Pose const & pose,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	scoring::ScoreFunctionOP const & scorefxn_sc, // Only hbond_sc_bb and hbond_sc
	scoring::ScoreFunctionOP const & scorefxn_bb, // Only hbond_lr_bb and hbond_sr_bb
	float atomic_depth_cutoff = 4.5f,
	float atomic_depth_probe_radius = 2.3f,
	float atomic_depth_resolution = 0.5f,
	float minimum_hb_cut = 0,
	bool all_atoms_active = false,
	float oversat_penalty = 1,
	bool assume_const_backbone = true
);

void
add_to_onebody(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	Size resnum,
	Size rotamer_id,
	float adder
);

struct ReweightData {

	ReweightData( pose::Pose const & pose_in, task::PackerTaskCOP const & task_in )
	: pose( pose_in ), task( task_in ), edge_reweights( task_in->IGEdgeReweights() )
	{}

	pose::Pose const & pose;
	task::PackerTaskCOP const & task;
	task::IGEdgeReweightContainerCOP edge_reweights;
	std::unordered_map< uint64_t, float > stored_edge_reweights;
};


void
add_to_twobody(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	ReweightData & reweight_data,
	Size resnum1,
	Size rotamer_id1,
	Size resnum2,
	Size rotamer_id2,
	float adder
);

uint64_t
map_key_twob( Size resnum1, Size resnum2 );

uint64_t
map_key_oneb( Size resnum1 );

}
}
}
}

#endif
