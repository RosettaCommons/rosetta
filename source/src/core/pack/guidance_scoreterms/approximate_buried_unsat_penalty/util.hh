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
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <basic/datacache/CacheableResRotPairFloatMap.hh>

#include <utility/vector1.hh>
#include <unordered_map>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {


struct UnsatCorrectionOptions {

public:
	UnsatCorrectionOptions() :
		nh2_wants_2(false),
		nh1_wants_1(false),
		hydroxyl_wants_h(false),
		carboxyl_wants_2(false)
	{}

	bool nh2_wants_2;
	bool nh1_wants_1;
	bool hydroxyl_wants_h;
	bool carboxyl_wants_2;
};

struct HBondBonusOptions {

public:
	HBondBonusOptions() :
		scorefxn_weight_(1),
		cross_chain_bonus_(0),
		ser_to_helix_bb_(0)
	{}

	float scorefxn_weight_;
	float cross_chain_bonus_;
	float ser_to_helix_bb_;

	bool any() const {
		return cross_chain_bonus_ != 0 || ser_to_helix_bb_ != 0;
	}

	float cross_chain_bonus() const {
		return cross_chain_bonus_ / scorefxn_weight_;
	}
	float ser_to_helix_bb() const {
		return ser_to_helix_bb_ / scorefxn_weight_;
	}
};


basic::datacache::CacheableResRotPairFloatMapOP
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
	bool assume_const_backbone = true,
	UnsatCorrectionOptions const & cor_opt = UnsatCorrectionOptions(),
	HBondBonusOptions const & bonus_opt = HBondBonusOptions()

);

void
add_to_onebody(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	utility::vector1<bool> const & is_asu,
	Size resnum,
	Size rotamer_id,
	float adder
);

struct ReweightData {

	ReweightData( pose::Pose const & pose_in, task::PackerTaskCOP const & task_in )
	: pose( pose_in ), task( task_in ), edge_reweights( task_in->IGEdgeReweights() )
	{}

	pose::Pose const & pose;
	task::PackerTaskCOP task;
	task::IGEdgeReweightContainerCOP edge_reweights;
	std::unordered_map< basic::datacache::ResRotPair, float, basic::datacache::ResRotPairHasher > stored_edge_reweights;
};


void
add_to_twobody(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	utility::vector1<bool> const & is_asu,
	ReweightData & reweight_data,
	Size resnum1,
	Size rotamer_id1,
	Size resnum2,
	Size rotamer_id2,
	float adder,
	core::conformation::symmetry::SymmetryInfoCOP const & symm_info /*can be nullptr*/
);


}
}
}
}

#endif
