// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.hh
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_scoring_lkball_LK_DomeHelper_hh
#define INCLUDED_core_scoring_lkball_LK_DomeHelper_hh

// Unit headers
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.fwd.hh>
#include <core/scoring/lkball/LK_DomeEnergy.hh>

// Package headers
#include <core/id/AtomID_Map.hh>
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


#define FULLY_OCCLUDED 127


// forward declaration for testing
class LK_DomeEnergyTests;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace lk_dome {



struct MyWaterHolder {
	bool is_bb;
	Real sol_value;
	Vector water_xyz;
	Vector base_xyz;
	Size iatom;

	MyWaterHolder( bool bb, Real sol, Vector const & wat, Vector const & base, Size iat=0 );

	bool operator<( MyWaterHolder const & o ) const;

	bool operator==( MyWaterHolder const & o ) const;
};

uint32_t
popcount_until_bit( uint32_t bit, uint32_t const * smallest_byte );

void
bit_set( int32_t bit, uint32_t * smallest_byte );

bool
bit_test( int32_t bit, uint32_t const * smallest_byte );


class LK_DomeHelper {
	friend class ::LK_DomeEnergyTests;

public:

	LK_DomeHelper(
		core::scoring::lkball::LK_DomeEnergyCOP const & lk_dome,
		Real lk_dome_weight,
		Real lk_dome_iso_weight,
		Real lk_dome_bridge_weight,
		Real lk_dome_bridge_uncpl_weight,
		Real lk_ball_bridge2_weight,
		Real lk_ball_bridge_uncpl2_weight
	);

	core::Real
	calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const &resvect,
		utility::vector1< core::Size > const & current_rotamer_ids,
		core::Size const substitution_position
	);

	void
	commit_considered_substitution();


	core::pack::rotamer_set::RotamerSets
	init_with_pose(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	Real current_score() const;



protected:

	void
	restore_from_shadow();

	void
	save_to_shadow();

	void
	reinit_with_resvect(
		utility::vector1< core::Size > const & current_rotamer_ids,
		bool skip_reset = false
	);

	void
	reset_calculation();

	void
	resize_arrays(
		core::pose::Pose const & pose
	);

	void
	init();

	void
	add_remove_rotamer(
		utility::vector1< core::Size > const & current_rotamer_ids,
		Size const substitution_position,
		bool add
	);

	void
	add_waters( Size position, Size rotamer );

	void
	init_water_scores( Size position );

	uint32_t
	get_score_data_offset( uint32_t water_seqpos, uint16_t water_id, uint32_t other_seqpos, uint32_t other_rotamer );

	void
	update_other_waters_inner( int8_t occlusion_delta, int16_t interact_delta, uint32_t water_seqpos, uint32_t total_water_offset );

	void
	update_other_waters( Size position, bool adding );

	void
	prepare_lkd_infos_and_assign_waters(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets,
		utility::vector1<utility::vector0<core::scoring::lkball::LKD_ResidueInfo>> & lkd_infos,
		utility::vector1<utility::vector0<MyWaterHolder>> & waters,
		utility::vector1<utility::vector0<utility::vector0<Size>>> & water_assignments
	);

	void
	prepare_simple_rotamer_values(
		utility::vector1<utility::vector0<MyWaterHolder>> & waters,
		utility::vector1<utility::vector0<utility::vector0<Size>>> const & water_assignments
	);

	void
	prepare_score_data(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets,
		utility::vector1<utility::vector0<core::scoring::lkball::LKD_ResidueInfo>> const & lkd_infos,
		utility::vector1<utility::vector0<MyWaterHolder>> const & waters,
		utility::vector1<utility::vector0<utility::vector0<Size>>> const & water_assignments
	);

	core::scoring::lkball::LK_DomeEnergyCOP lk_dome() const;

	Size
	memory_use() const;


private:

	core::scoring::ScoreFunctionOP fake_lr_scorefxn_;
	core::scoring::lkball::LK_DomeEnergyCOP lk_dome_;

	Real lk_dome_weight_;
	Real lk_dome_iso_weight_;
	Real lk_dome_bridge_weight_;
	Real lk_dome_bridge_uncpl_weight_;
	Real lk_ball_bridge2_weight_;
	Real lk_ball_bridge_uncpl2_weight_;


private:

	// runtime data

	int64_t total_score_;

	// utility::vector1<bool> score_needs_recalc_;

	utility::vector1<Size> current_rotamer_at_seqpos_;
	utility::vector1<int32_t> current_score_at_seqpos_;

	// per-seqpos data
	// should be vector of vectors but this is better cache-performance
	// size = pose.size() * waters_per_seqpos_
	utility::vector1<int16_t> current_waters_;
	utility::vector1<int16_t> current_water_occl_sum_;
	utility::vector1<int32_t> current_water_interact_sum_;


	utility::vector1<bool> shadow_mismatch_;

	int64_t shadow_total_score_;
	utility::vector1<Size> shadow_current_rotamer_at_seqpos_;
	utility::vector1<int32_t> shadow_current_score_at_seqpos_;
	utility::vector1<int16_t> shadow_current_waters_;
	utility::vector1<int16_t> shadow_current_water_occl_sum_;
	utility::vector1<int32_t> shadow_current_water_interact_sum_;



	// permanent data

	Real partial_score_scaler_;
	Real score_scaler_;
	uint8_t occlusion_min_;
	uint8_t occlusion_span_;

	// used for indexing current_* arrays
	uint32_t waters_per_seqpos_;

	// Lookups into other arrays
	utility::vector1<uint32_t> seqpos_to_rotamer_to_waters_masks_;
	utility::vector1<uint32_t> seqpos_to_water_to_pwp_index_list_;
	utility::vector1<uint32_t> seqpos_to_riww_at_pos_index_list_;

	// How many unique waters at this position
	utility::vector1<uint16_t> seqpos_num_unique_waters_;

	// Uber hot
	// number of rotamers divided by 32 and rounded up
	utility::vector1<uint8_t> seqpos_num_rotamers_32_;

	// stores masks of which waters are active for this rotamer
	//   -- masks are different for each seqpos
	// major indexed by: seqpos_to_rotamer_to_waters_masks_
	// minor indexed by: irot * seqpos_num_unique_waters_ // 32
	// data size: seqpos_num_unique_waters_ // 32
	// data: masks of which waters to use
	// -- least significant uint32 first
	utility::vector1<uint32_t> rotamer_to_waters_masks_;


	// stores index and length of data in water_to_pwp_index_list
	// major indexed by: seqpos_to_water_to_pwp_index_list_
	// minor indexed by: water_id * 2
	// data size: 2
	// data: idx into pos_water_pair_to_score
	//       num elements in pos_water_pair_to_score
	utility::vector1<uint32_t> water_to_pwp_index_list_;

	#define PWP_OTHER_SEQPOS 0
	//#define PWP_number_of_bitfields 1
	#define PWP_OFFSET_INTO_SCORE_DATA 1
	#define PWP_BITFIELDS 2

	// stores all rotamer interactions at seqpos with a specific water
	// major indexed by: water_to_pwp_index_list_
	// data size: num elements given by water_to_pwp_index_list_
	//               each element is variable size
	// data: uint32_t other_seqpos
	//       uint32_t offset_into_score_data
	//       uint32_t[] bitfields -- size of bitfields defined by seqpos_num_rotamers_32_
	// -- least significant uint32_t first
	utility::vector1<uint32_t> pos_water_pair_to_score_;


	// Stores occlusion fraction and interaction score each encoded in 8 bits
	// major indexed by: pos_water_pair_to_score_ -> offset_into_score_data
	// minor indexed by: popcount(less-significant bitfields) * 3
	// data size: 3
	// data: int8 occlusion_fraction
	//       int16 interaction score
	utility::vector1<int8_t> score_data_;


	// Stores index and length for a rotamer into riww_at_pos
	// major indexed by: seqpos_to_riww_at_pos_index_list_
	// minor indexed by: irot * 2
	// data size: 2
	// data: idx into riww_at_pos
	//       size of entry in riww_at_pos
	utility::vector1<uint32_t> riww_at_pos_index_list_;


	// Stores list of seqpos and water masks for each rotamer.
	//    Used to see who this rotamer occludes and scores with
	// major indexed by: riww_at_pos_index_list_
	// data size: elements given by riww_at_pos_index_list_. each element variable
	// data: uint32_t water_seqpos
	//       uint32_t[] masks of water interactions -- length defined by seqpos_num_unique_waters_
	utility::vector1<uint32_t> rotamer_interacts_with_waters_at_pos_;



	// debugging stuff

private:

	bool debug_store_waters_;
	utility::pointer::shared_ptr< utility::vector1<utility::vector0<MyWaterHolder>> > debug_waters_;


};

}
} //guidance_scoreterms
} //pack
} //core

#endif
