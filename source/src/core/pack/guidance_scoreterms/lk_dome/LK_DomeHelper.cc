// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.cc
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.hh>


// Package headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/conformation/util.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/select/util.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>

// Basic Headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <utility/graph/Graph.hh>

#include <bitset>

static basic::Tracer TR("core.scoring.lkball.LK_DomeHelper");

using namespace core::scoring;
using namespace core::scoring::lkball;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace lk_dome {





LK_DomeHelper::LK_DomeHelper(
	core::scoring::lkball::LK_DomeEnergyCOP const & lk_dome,
	Real lk_dome_weight,
	Real lk_dome_iso_weight,
	Real lk_dome_bridge_weight,
	Real lk_dome_bridge_uncpl_weight,
	Real lk_ball_bridge2_weight,
	Real lk_ball_bridge_uncpl2_weight
) :
	lk_dome_( lk_dome ),
	lk_dome_weight_( lk_dome_weight ),
	lk_dome_iso_weight_( lk_dome_iso_weight ),
	lk_dome_bridge_weight_( lk_dome_bridge_weight ),
	lk_dome_bridge_uncpl_weight_( lk_dome_bridge_uncpl_weight ),
	lk_ball_bridge2_weight_( lk_ball_bridge2_weight ),
	lk_ball_bridge_uncpl2_weight_( lk_ball_bridge_uncpl2_weight ),
	debug_store_waters_( false )
{
	init();
}


bool
bit_test( int32_t bit, uint32_t const * smallest_byte ) {
	int32_t offset = bit >> 5;
	int32_t mask = bit & 0x1F;
	return (*(smallest_byte + offset)) & ( 1 << mask );
}

void
bit_set( int32_t bit, uint32_t * smallest_byte ) {
	int32_t offset = bit >> 5;
	int32_t mask = bit & 0x1F;
	(*(smallest_byte + offset)) |= ( 1 << mask );
}


uint32_t
popcount_until_bit( uint32_t bit, uint32_t const * smallest_byte ) {
	uint32_t final_mask = uint64_t(0xFFFFFFFF) >> ( 32 - (bit & 0x1F) ); // can't shift a uint32_t by 32, I tried
	uint32_t whole_bytes = bit >> 5;

	uint32_t sum = 0;
	while ( whole_bytes > 0 ) {
		sum += std::bitset<32>( *smallest_byte ).count();
		smallest_byte++;
		whole_bytes--;
	}
	sum += std::bitset<32>( (*smallest_byte) & final_mask ).count();
	return sum;
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// PACKER RUNTIME ROUTINES /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// This is the main function that gets called during packing. Notably this function either
//  adjusts the current internal resvector or restarts from scratch.
//  Before every substition, the shadow is copied back into our working registers so that
//  we can easily recover from failed packer moves.
core::Real
LK_DomeHelper::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &,
	utility::vector1< core::Size > const & current_rotamer_ids,
	core::Size const substitution_position
) {

	if ( substitution_position == 0 ) {
		reinit_with_resvect( current_rotamer_ids );
	} else {
		restore_from_shadow();

		add_remove_rotamer( current_rotamer_at_seqpos_, substitution_position, false );
		add_remove_rotamer( current_rotamer_ids, substitution_position, true );


	}

	return total_score_ * score_scaler_;
}

// All we need to do to commit a substitution is to save our current registers into the shadow
void
LK_DomeHelper::commit_considered_substitution() {
	// std::cout << "   commit" << std::endl;
	save_to_shadow();
}

// This function copies all of the shadow registers back into the working registers. This is used
//  so that we can ignore a failed packer substitution
void
LK_DomeHelper::restore_from_shadow() {

	for ( Size seqpos = 1; seqpos <= current_rotamer_at_seqpos_.size(); seqpos ++ ) {
		if ( ! shadow_mismatch_[seqpos] ) continue;
		shadow_mismatch_[seqpos] = false;

		current_rotamer_at_seqpos_[seqpos] = shadow_current_rotamer_at_seqpos_[seqpos];
		current_score_at_seqpos_[seqpos] = shadow_current_score_at_seqpos_[seqpos];

		Size water_offset = (seqpos - 1)*waters_per_seqpos_ + 1;
		for ( Size iwat = 0; iwat < waters_per_seqpos_; iwat++ ) {
			current_waters_[ water_offset + iwat ] = shadow_current_waters_[ water_offset + iwat ];
			current_water_occl_sum_[ water_offset + iwat ] = shadow_current_water_occl_sum_[ water_offset + iwat ];
			current_water_interact_sum_[ water_offset + iwat ] = shadow_current_water_interact_sum_[ water_offset + iwat ];
		}

	}

	total_score_ = shadow_total_score_;


}

// Save our corrent working registers into the shadow. This is what it means to "commit" a change.
void
LK_DomeHelper::save_to_shadow() {

	for ( Size seqpos = 1; seqpos <= current_rotamer_at_seqpos_.size(); seqpos ++ ) {
		if ( ! shadow_mismatch_[seqpos] ) continue;
		shadow_mismatch_[seqpos] = false;

		shadow_current_rotamer_at_seqpos_[seqpos] = current_rotamer_at_seqpos_[seqpos];
		shadow_current_score_at_seqpos_[seqpos] = current_score_at_seqpos_[seqpos];

		Size water_offset = (seqpos - 1)*waters_per_seqpos_ + 1;
		for ( Size iwat = 0; iwat < waters_per_seqpos_; iwat++ ) {
			shadow_current_waters_[ water_offset + iwat ] = current_waters_[ water_offset + iwat ];
			shadow_current_water_occl_sum_[ water_offset + iwat ] = current_water_occl_sum_[ water_offset + iwat ];
			shadow_current_water_interact_sum_[ water_offset + iwat ] = current_water_interact_sum_[ water_offset + iwat ];
		}

	}

	shadow_total_score_ = total_score_;

}


void
LK_DomeHelper::add_remove_rotamer(
	utility::vector1< core::Size > const & current_rotamer_ids,
	Size const substitution_position,
	bool add
) {
	shadow_mismatch_[ substitution_position ] = true;

	if ( add ) {
		current_rotamer_at_seqpos_[ substitution_position ] = current_rotamer_ids[ substitution_position ];

		// first we clear out the old water scores because we're going to nuke them
		total_score_ -= current_score_at_seqpos_[ substitution_position ];
		current_score_at_seqpos_[ substitution_position ] = 0;

		// next, we identify the new waters and add their scores
		add_waters( substitution_position, current_rotamer_ids[ substitution_position ] );
		init_water_scores( substitution_position );

		// lastly, we update waters at other positions that we are occluding
		update_other_waters( substitution_position, true );


	} else {

		// In remove, we're primarily concerned with removing our occlusions. This is the only thing rotamer specific
		update_other_waters( substitution_position, false );
	}

}


// All the int16_t scores start out invalid when this is called
void
LK_DomeHelper::add_waters( Size position, Size rotamer ) {

	uint32_t water_offset = (position-1)*waters_per_seqpos_ + 1;

	uint16_t this_position_unique_waters = seqpos_num_unique_waters_[ position ];
	if ( this_position_unique_waters == 0 ) {
		for ( Size i = 0; i < waters_per_seqpos_; i++ ) {
			current_waters_[water_offset + i] = -1;
		}
		return;
	}

	uint16_t this_position_unique_waters_32 = ( (this_position_unique_waters-1) >> 5) + 1;

	uint32_t major_index = seqpos_to_rotamer_to_waters_masks_[ position ];
	uint32_t minor_index = rotamer * this_position_unique_waters_32;
	uint32_t mask_cursor = major_index + minor_index;

	Size iwat = 0;
	Size bit = 0;
	for ( Size ibyte = 0; ibyte < this_position_unique_waters_32; ibyte++ ) {
		uint32_t test_mask = 1;
		uint32_t mask = rotamer_to_waters_masks_[ mask_cursor + ibyte ];
		for ( Size inner_bit = 0; inner_bit < 32; inner_bit++ ) {
			if ( (mask & test_mask) != 0 ) {
				current_waters_[ water_offset + iwat ] = bit;
				iwat++;
			}
			bit++;
			test_mask <<= 1;
			if ( bit >= this_position_unique_waters ) break;
		}
	}
	Size num_waters = iwat;
	debug_assert( num_waters <= waters_per_seqpos_ );
	// Clean out the ones we're not going to get to
	for ( iwat = num_waters; iwat < waters_per_seqpos_; iwat++ ) {
		current_waters_[water_offset + iwat] = -1;
	}
}

void
LK_DomeHelper::init_water_scores( Size position ) {

	shadow_mismatch_[ position ] = true;

	uint32_t water_offset = (position-1)*waters_per_seqpos_ + 1;

	uint32_t major_index = seqpos_to_water_to_pwp_index_list_[ position ];

	// Loop over waters at this position
	for ( Size iwat = 0; iwat < waters_per_seqpos_; iwat++ ) {
		int16_t water_id = current_waters_[ water_offset + iwat ];
		if ( water_id == -1 ) {
			current_water_occl_sum_[ water_offset + iwat ] = 0;
			current_water_interact_sum_[ water_offset + iwat ] = 0;
			continue;
		}

		uint32_t minor_index = water_id * 2;

		uint32_t pwp_cursor = water_to_pwp_index_list_[major_index + minor_index];
		uint32_t pwp_elements = water_to_pwp_index_list_[major_index + minor_index + 1];

		int16_t this_water_occl_sum = 0;
		int32_t this_water_interact_sum = 0;

		// Loop over pwp looking to see if the current rotamer interacts with this water
		for ( Size i_element = 1; i_element <= pwp_elements; i_element++ ) {

			uint32_t other_seqpos = pos_water_pair_to_score_[ pwp_cursor + PWP_OTHER_SEQPOS ];
			Size other_rotamer = current_rotamer_at_seqpos_[ other_seqpos ];

			// If set, the current rotamer interacts with this water
			if ( bit_test( (int32_t)other_rotamer, &pos_water_pair_to_score_[ pwp_cursor + PWP_BITFIELDS ] ) ) {

				uint32_t offset_into_score_data = pos_water_pair_to_score_[ pwp_cursor + PWP_OFFSET_INTO_SCORE_DATA ];
				uint32_t prev_popcount = popcount_until_bit( other_rotamer, &pos_water_pair_to_score_[ pwp_cursor + PWP_BITFIELDS ] );

				// * 3 because each entry in score_data_ is 3 bytes
				uint32_t final_offset_into_score_data = offset_into_score_data + prev_popcount * 3;

				// We get both the lk_dome scores of that rotamer to our water as well as the other rotamer's
				//   occlusion on this water
				int8_t occlusion = score_data_[ final_offset_into_score_data ];
				int16_t score_val = *((int16_t *)&score_data_[ final_offset_into_score_data + 1 ]);

				this_water_occl_sum += occlusion;
				this_water_interact_sum += score_val;
			}

			pwp_cursor += PWP_BITFIELDS + seqpos_num_rotamers_32_[ other_seqpos ];
		}

		current_water_occl_sum_[ water_offset + iwat ] = this_water_occl_sum;
		current_water_interact_sum_[ water_offset + iwat ] = this_water_interact_sum;

		// Ok, now we have accounted for all other rotamers acting on us. Tally up the scores

		// Fully occluded, don't update scores
		if ( this_water_occl_sum >= FULLY_OCCLUDED ) continue;

		uint32_t avail = occlusion_span_;
		if ( this_water_occl_sum > occlusion_min_ ) {
			avail = occlusion_span_ - (this_water_occl_sum - occlusion_min_);

			debug_assert( occlusion_span_ >= (this_water_occl_sum - occlusion_min_) );
		}

		// uint32_t avail = FULLY_OCCLUDED - this_water_occl_sum;
		int32_t delta_score = avail * this_water_interact_sum;

		total_score_ += delta_score;
		current_score_at_seqpos_[ position ] += delta_score;
	}
}



// Do the crazy lookup into score_data_
uint32_t
LK_DomeHelper::get_score_data_offset( uint32_t water_seqpos, uint16_t water_id, uint32_t other_seqpos, uint32_t other_rotamer ) {
	uint32_t major_offset = seqpos_to_water_to_pwp_index_list_[ water_seqpos ];
	uint32_t minor_offset = water_id * 2;

	uint32_t pwp_cursor = water_to_pwp_index_list_[ major_offset + minor_offset ];
	uint32_t pwp_elements = water_to_pwp_index_list_[ major_offset + minor_offset + 1];

	// surf through pwp until we hit our seqpos
	// This is potentially slow, but idk
	for ( uint32_t i_element = 1; i_element <= pwp_elements; i_element++ ) {
		uint32_t pwp_other_seqpos = pos_water_pair_to_score_[ pwp_cursor + PWP_OTHER_SEQPOS ];
		if ( pwp_other_seqpos != other_seqpos ) {
			uint32_t num_bitfields = seqpos_num_rotamers_32_[ pwp_other_seqpos ];
			pwp_cursor += PWP_BITFIELDS + num_bitfields;
			continue;
		}

		// If you got this far, you should probably be looking up a valid score...
		debug_assert( bit_test( other_rotamer, &pos_water_pair_to_score_[ pwp_cursor + PWP_BITFIELDS ] ) );

		uint32_t offset_into_score_data = pos_water_pair_to_score_[ pwp_cursor + PWP_OFFSET_INTO_SCORE_DATA ];
		uint32_t prev_popcount = popcount_until_bit( other_rotamer, &pos_water_pair_to_score_[ pwp_cursor + PWP_BITFIELDS ] );

		// Each entry in score_data is 3 bytes so multiply by 3
		return offset_into_score_data + prev_popcount * 3;

	}

	utility_exit_with_message("LK_DomeHelper: Error in get_score_data_offset");
}

// Update all scores as a result of an occlusion changing
// This is sort of complicated because the occlusion clips to FULLY_OCCLUDED
// Also, don't overflow your integers...
void
LK_DomeHelper::update_other_waters_inner( int8_t occlusion_delta, int16_t interact_delta, uint32_t water_seqpos, uint32_t total_water_offset ) {
	int16_t old_occlusion = current_water_occl_sum_[total_water_offset];
	int16_t new_occlusion = old_occlusion + (int16_t)occlusion_delta;
	current_water_occl_sum_[total_water_offset] = new_occlusion;

	shadow_mismatch_[ water_seqpos ] = true;

	int32_t real_delta = occlusion_delta;
	if ( occlusion_delta > 0 ) {
		if ( old_occlusion >= FULLY_OCCLUDED ) {
			current_water_interact_sum_[ total_water_offset ] += interact_delta;
			return; // scores don't change because we were already maxed
		}
		if ( new_occlusion <= occlusion_min_ ) {
			// We're still maxing avail. However, this means scores do change
			real_delta = 0;
		} else {
			if ( old_occlusion <= occlusion_min_ ) {
				if ( new_occlusion >= FULLY_OCCLUDED ) {
					real_delta = occlusion_span_; // max possible movement
				} else {
					real_delta = new_occlusion - occlusion_min_;
				}
			} else {
				if ( new_occlusion >= FULLY_OCCLUDED ) {
					real_delta = FULLY_OCCLUDED - old_occlusion; // just move to max
				}
			}
		}
	} else {
		if ( new_occlusion >= FULLY_OCCLUDED ) {
			current_water_interact_sum_[ total_water_offset ] += interact_delta;
			return; // scores don't change because we're still maxed
		}
		if ( old_occlusion <= occlusion_min_ ) {
			// We're still maxing avail. However, this means scores do change
			real_delta = 0;
		} else {
			if ( new_occlusion <= occlusion_min_ ) {
				// Crossing the min border
				if ( old_occlusion >= FULLY_OCCLUDED ) {
					// completely crossing the linear range
					real_delta = -occlusion_span_;
				} else {
					// started in linear and are now below
					real_delta = occlusion_min_ - old_occlusion;
				}
			} else {
				// We're not touching the min border
				if ( old_occlusion >= FULLY_OCCLUDED ) {
					// previously we were above max
					real_delta = new_occlusion - FULLY_OCCLUDED;
				}
			}
		}

	}
	// Negative sign here because if occlusion goes up, score goes down
	int32_t delta_score = -real_delta * (int32_t)current_water_interact_sum_[ total_water_offset ];

	// Ok, now we change the interaction value
	current_water_interact_sum_[ total_water_offset ] += interact_delta;
	if ( new_occlusion < FULLY_OCCLUDED ) {

		uint32_t avail = occlusion_span_;
		if ( new_occlusion > occlusion_min_ ) {
			avail = occlusion_span_ - (new_occlusion - occlusion_min_);

			debug_assert( occlusion_span_ >= (new_occlusion - occlusion_min_) );
		}

		delta_score += avail * interact_delta;
	}

	current_score_at_seqpos_[ water_seqpos ] += delta_score;
	total_score_ += delta_score;
}


void
LK_DomeHelper::update_other_waters( Size position, bool adding ) {

	uint32_t major_index = seqpos_to_riww_at_pos_index_list_[ position ];
	uint32_t minor_index = current_rotamer_at_seqpos_[ position ] * 2;

	uint32_t cursor = riww_at_pos_index_list_[ major_index + minor_index ];
	uint32_t num_positions   = riww_at_pos_index_list_[ major_index + minor_index + 1 ];

	for ( Size ipos = 0; ipos < num_positions; ipos++ ) {
		uint32_t water_seqpos = rotamer_interacts_with_waters_at_pos_[ cursor ];
		uint32_t * bitmask = &rotamer_interacts_with_waters_at_pos_[ cursor + 1 ];

		uint32_t water_offset = (water_seqpos-1)*waters_per_seqpos_ + 1;
		for ( Size iwat = 0; iwat < waters_per_seqpos_; iwat++ ) {
			int16_t wat_id = current_waters_[water_offset + iwat];
			if ( wat_id == -1 ) continue;
			if ( ! bit_test( wat_id, bitmask ) ) continue;

			// Ok, now we know that this rotamer occludes this water to some degree
			uint32_t score_data_offset = get_score_data_offset( water_seqpos, wat_id, position, current_rotamer_at_seqpos_[ position ] );

			int8_t occlusion = score_data_[score_data_offset];
			int16_t interact = *((int16_t *)&score_data_[score_data_offset + 1]);

			if ( adding ) {
				update_other_waters_inner( occlusion, interact, water_seqpos, water_offset + iwat );
			} else {
				update_other_waters_inner(-occlusion, -interact, water_seqpos, water_offset + iwat );
			}
		}
		cursor += 1 + ( ( seqpos_num_unique_waters_[water_seqpos] - 1 ) >> 5 ) + 1;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// END PACKER RUNTIME ROUTINES ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


// This function resets the calculation with a fresh set of rotamers. We simply add each residue
//  one by one and then recalculate sasa and saps at the end
void
LK_DomeHelper::reinit_with_resvect(
	utility::vector1< core::Size > const & current_rotamer_ids,
	bool skip_reset /* = false */
) {
	if ( ! skip_reset ) reset_calculation();

	// Ok, we do this in 2 stages

	// Stage 1. Add waters everywhere but don't touch the scores
	for ( Size seqpos = 1; seqpos <= current_rotamer_ids.size(); seqpos++ ) {
		current_rotamer_at_seqpos_[seqpos] = current_rotamer_ids[seqpos];
		add_waters( seqpos, current_rotamer_ids[ seqpos ] );
	}

	// Stage 2. Have each water accumulate scores. They each take care of their own
	//  occlusion and don't touch anyone else

	for ( Size seqpos = 1; seqpos <= current_rotamer_ids.size(); seqpos++ ) {
		init_water_scores( seqpos ); // sets shadow_mismatch true
	}

	// no need to call update_other_waters() because the waters handled it themselves

	save_to_shadow();
}


// Clear and resize all arrays that accumulate things related to scores.
// This will not prepare for a new pose however
void
LK_DomeHelper::reset_calculation() {

	Size pose_size = seqpos_to_rotamer_to_waters_masks_.size();

	total_score_ = 0;

	current_rotamer_at_seqpos_.clear();
	current_rotamer_at_seqpos_.resize( pose_size );

	current_score_at_seqpos_.clear();
	current_score_at_seqpos_.resize( pose_size );

	current_waters_.clear();
	current_waters_.resize( pose_size * waters_per_seqpos_, -1 );
	current_water_occl_sum_.clear();
	current_water_occl_sum_.resize( pose_size * waters_per_seqpos_, 0 );
	current_water_interact_sum_.clear();
	current_water_interact_sum_.resize( pose_size * waters_per_seqpos_, 0 );


	shadow_mismatch_.clear();
	shadow_mismatch_.resize( pose_size );

	shadow_total_score_ = 0;

	shadow_current_rotamer_at_seqpos_.clear();
	shadow_current_rotamer_at_seqpos_.resize( pose_size );

	shadow_current_score_at_seqpos_.clear();
	shadow_current_score_at_seqpos_.resize( pose_size );

	shadow_current_waters_.clear();
	shadow_current_waters_.resize( pose_size * waters_per_seqpos_, -1 );
	shadow_current_water_occl_sum_.clear();
	shadow_current_water_occl_sum_.resize( pose_size * waters_per_seqpos_, 0 );
	shadow_current_water_interact_sum_.clear();
	shadow_current_water_interact_sum_.resize( pose_size * waters_per_seqpos_, 0 );
}



// Resize all arrays to the size of a new pose
void
LK_DomeHelper::resize_arrays(
	core::pose::Pose const &
) {



}


// Called immediately after construction
// This just sets up things that don't change based on input pose
void
LK_DomeHelper::init() {


	fake_lr_scorefxn_ = utility::pointer::make_shared<scoring::ScoreFunction>();
	fake_lr_scorefxn_->set_weight(scoring::fa_rep, 1);
	scoring::methods::EnergyMethodOptions opts = fake_lr_scorefxn_->energy_method_options();
	opts.etable_options().max_dis = lk_dome_->atomic_interaction_cutoff();
	fake_lr_scorefxn_->set_energy_method_options( opts );

}


MyWaterHolder::MyWaterHolder( bool bb, Real sol, Vector const & wat, Vector const & base, Size iat )
: is_bb( bb ),
	sol_value( sol),
	water_xyz( wat ),
	base_xyz( base ),
	iatom( iat )
{}

bool
MyWaterHolder::operator<( MyWaterHolder const & o ) const {
	if ( is_bb != o.is_bb ) return is_bb > o.is_bb; // we want bb atoms to come first for water-matching during packing

	if ( sol_value != o.sol_value ) return sol_value < o.sol_value;

	// idk about the float thing. Ask frank
	for ( int i = 0; i < 3; i++ ) {
		if ( float( water_xyz[i]) != float(o.water_xyz[i]) ) {
			return water_xyz[i] < o.water_xyz[i];
		}
	}

	// idk about the float thing. Ask frank
	for ( int i = 0; i < 3; i++ ) {
		if ( float( base_xyz[i]) != float(o.base_xyz[i]) ) {
			return base_xyz[i] < o.base_xyz[i];
		}
	}

	return false;
}

bool
MyWaterHolder::operator==( MyWaterHolder const & o ) const {
	if ( is_bb != o.is_bb ) return false;

	if ( sol_value != o.sol_value ) return false;

	// idk about the float thing. Ask frank
	for ( int i = 0; i < 3; i++ ) {
		if ( float( water_xyz[i]) != float(o.water_xyz[i]) ) {
			return false;
		}
	}

	// idk about the float thing. Ask frank
	for ( int i = 0; i < 3; i++ ) {
		if ( float( base_xyz[i]) != float(o.base_xyz[i]) ) {
			return false;
		}
	}

	return true;
}


void
LK_DomeHelper::prepare_lkd_infos_and_assign_waters(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1<utility::vector0<LKD_ResidueInfo>> & lkd_infos,
	utility::vector1<utility::vector0<MyWaterHolder>> & waters,
	utility::vector1<utility::vector0<utility::vector0<Size>>> & water_assignments
) {
	lkd_infos.resize( pose.size() );
	waters.resize( pose.size() );
	water_assignments.resize( pose.size() );

	pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

		pack::rotamer_set::RotamerSetCOP const & rotset = rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
			rotamer_sets.rotamer_set_for_residue( seqpos )
			: fake_rotset;

		lkd_infos[seqpos].reserve( rotset->num_rotamers() + 1 );

		// std::map so that we can avoid duplicating waters
		std::map< MyWaterHolder, Size > found_waters;

		// Store this so we don't have to iterate the entire pose in loop 3
		utility::vector0<utility::vector0<MyWaterHolder>> waters_on_res( rotset->num_rotamers() + 1);

		for ( Size irot = 0; irot <= rotset->num_rotamers(); irot++ ) {
			conformation::Residue const & res = irot == 0 ? pose.residue(seqpos) : *rotset->rotamer( irot );
			lkd_infos[seqpos].emplace_back( res, &*lk_dome_ );

			LKD_ResidueInfo const & info = lkd_infos[seqpos][irot];

			if ( ! info.has_waters() ) continue;

			utility::vector1< Real > const & sol_values = info.water_sol_values();
			WaterCoords const & water_coords = info.waters();
			utility::vector1< Size > const & n_attached_waters = info.n_attached_waters();


			for ( Size iatom = 1; iatom <= res.nheavyatoms(); iatom++ ) {
				Size atom_offset = info.water_offset_for_atom()[iatom];
				for ( Size iwat = 1; iwat <= n_attached_waters[iatom]; iwat++ ) {
					found_waters.emplace(std::piecewise_construct,
						std::forward_as_tuple( res.atom_is_backbone(iatom), sol_values[iatom], water_coords[iwat + atom_offset], res.xyz(iatom), iatom),
						std::forward_as_tuple(0)
					);
					waters_on_res[irot].emplace_back( res.atom_is_backbone(iatom), sol_values[iatom], water_coords[iwat + atom_offset], res.xyz(iatom), iatom );
				}
			}
		}

		// So unfortunately, floats can't be hashed and std::set/std::map sorts the elements as you add them
		// We have to use 3 for-loops to do the assignment instead of 1

		for ( auto & pair : found_waters ) {
			pair.second = waters[seqpos].size();
			waters[seqpos].push_back( pair.first );
		}

		// This can easily be doubled by using uint16_t.max_value() instead of -1 to denote no water
		runtime_assert( waters[seqpos].size() < 32768 );

		water_assignments[seqpos].resize( rotset->num_rotamers() + 1 );

		// The waters have been numbered, assign them back
		for ( Size irot = 0; irot < waters_on_res.size(); irot++ ) {
			for ( Size iwat = 0; iwat < waters_on_res[irot].size(); iwat++ ) {
				water_assignments[seqpos][irot].push_back( found_waters.at( waters_on_res[irot][iwat] ) );
			}
		}
	}
}


void
LK_DomeHelper::prepare_simple_rotamer_values(
	utility::vector1<utility::vector0<MyWaterHolder>> & waters,
	utility::vector1<utility::vector0<utility::vector0<Size>>> const & water_assignments
) {
	seqpos_to_rotamer_to_waters_masks_.clear();
	seqpos_num_unique_waters_.clear();
	seqpos_num_rotamers_32_.clear();
	rotamer_to_waters_masks_.clear();

	Size max_waters_on_rotamer = 0;

	for ( Size seqpos = 1; seqpos <= waters.size(); seqpos++ ) {


		// Divide by 32
		Size num_rots = water_assignments[seqpos].size();
		runtime_assert( num_rots < 255*32 ); // 8160 // could be increased without much issue
		seqpos_num_rotamers_32_.push_back( uint8_t( ( (num_rots-1) >> 5) + 1 ) );

		// Simple progressive index thing
		seqpos_to_rotamer_to_waters_masks_.push_back( rotamer_to_waters_masks_.size() + 1 );

		Size unique_waters = waters[seqpos].size();
		seqpos_num_unique_waters_.push_back( unique_waters );
		if ( unique_waters == 0 ) continue;
		Size unique_waters_32 = ( (unique_waters-1) >> 5) + 1;

		for ( Size irot = 0; irot < num_rots; irot++ ) {

			rotamer_to_waters_masks_.resize( rotamer_to_waters_masks_.size() + unique_waters_32, 0 );

			// lol, watch out, the resize can reallocate the vector
			uint32_t * smallest_byte = &rotamer_to_waters_masks_[ rotamer_to_waters_masks_.size() - (unique_waters_32 - 1)];

			Size waters_on_rotamer = water_assignments[seqpos][irot].size();

			// Set the bits denoting our waters
			for ( Size iwat = 0; iwat < waters_on_rotamer; iwat++ ) {
				Size water_id = water_assignments[seqpos][irot][iwat];

				debug_assert( water_id < unique_waters_32 * 32 ); // the next function will corrupt memory without checking
				bit_set( uint32_t(water_id), smallest_byte );
			}

			max_waters_on_rotamer = std::max<Size>( max_waters_on_rotamer, waters_on_rotamer );
		}
	}

	waters_per_seqpos_ = max_waters_on_rotamer;
}


void
LK_DomeHelper::prepare_score_data(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1<utility::vector0<LKD_ResidueInfo>> const & lkd_infos,
	utility::vector1<utility::vector0<MyWaterHolder>> const & waters,
	utility::vector1<utility::vector0<utility::vector0<Size>>> const & water_assignments
) {
	// This is the quickest way to get a graph of residues within the interaction distance
	pack::task::PackerTaskOP task = rotamer_sets.task()->clone();
	task->clean_residue_task( pose.residue(1), 1, pose );
	utility::vector1< Distance > residue_radii = pack::find_residue_max_radii( pose, task );
	utility::graph::GraphOP graph = create_packer_graph( pose, *fake_lr_scorefxn_, task, pose.size(), residue_radii );

	Real longest_radius = *std::max_element( residue_radii.begin(), residue_radii.end() );
	Real water_atom_range = lk_dome_->water_atom_interaction_cutoff();
	Real water_atom_range2 = water_atom_range * water_atom_range;
	Real longest_cb_water_interaction = longest_radius + water_atom_range;
	Real longest_cb_water_interaction2 = longest_cb_water_interaction * longest_cb_water_interaction;


	pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();


	utility::vector1<Real> real_score_data;
	seqpos_to_water_to_pwp_index_list_.clear();
	water_to_pwp_index_list_.clear();
	pos_water_pair_to_score_.clear();


	Real largest_single_interaction_mag = 0;

	utility::vector1<Size> num_waters_32s(pose.size());
	for ( Size water_seqpos = 1; water_seqpos <= pose.size(); water_seqpos++ ) {
		Size num_waters = waters[water_seqpos].size();
		if ( num_waters == 0 ) {
			num_waters_32s[water_seqpos] = 0;
		} else {
			num_waters_32s[water_seqpos] = ( (num_waters-1) >> 5) + 1;
		}
	}

	// This is the underlying data for riww
	utility::vector1< utility::vector0< utility::vector1< utility::vector0<uint32_t>>>> oseqpos_orotamer_wseqpos_water_interact_mask( pose.size() );
	for ( Size other_seqpos = 1; other_seqpos <= pose.size(); other_seqpos++ ) {

		Size irots_at_other = water_assignments[other_seqpos].size();

		oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos].resize(irots_at_other);
		for ( Size irot = 0; irot < irots_at_other; irot++ ) {

			oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos][irot].resize( pose.size() );

			for ( Size water_seqpos = 1; water_seqpos <= pose.size(); water_seqpos++ ) {
				oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos][irot][water_seqpos].resize( num_waters_32s[water_seqpos] );
			}
		}
	}

	// Here we get all interactions between all waters and all rotamers
	// We're filling in pwp as this moves along
	// score_data needs to be turned into int8_t after this function
	for ( Size water_seqpos = 1; water_seqpos <= pose.size(); water_seqpos++ ) {

		seqpos_to_water_to_pwp_index_list_.push_back( water_to_pwp_index_list_.size() + 1 );

		Size num_waters = waters[water_seqpos].size();
		if ( num_waters == 0 ) continue;


		for ( Size water_id = 0; water_id < num_waters; water_id++ ) {

			// push back index
			water_to_pwp_index_list_.push_back(pos_water_pair_to_score_.size() + 1);

			MyWaterHolder const & wat = waters[water_seqpos][water_id];

			Size num_seqpos_added = 0;

			for ( auto
					iter = graph->get_node( water_seqpos )->const_edge_list_begin(),
					iter_end = graph->get_node( water_seqpos )->const_edge_list_end();
					iter != iter_end; ++iter
					) {
				Size other_seqpos = (*iter)->get_other_ind(water_seqpos);
				Size seqpos_num_rotamers_32 = seqpos_num_rotamers_32_[ other_seqpos ];

				bool is_protein = pose.residue(other_seqpos).is_protein();
				Real nbr_dist2 = pose.residue(other_seqpos).nbr_atom_xyz().distance_squared( wat.water_xyz );
				if ( nbr_dist2 > longest_cb_water_interaction2 ) continue;

				utility::vector0<uint32_t> this_mask(seqpos_num_rotamers_32);
				bool any_interactions = false;
				Size starting_size = real_score_data.size();

				pack::rotamer_set::RotamerSetCOP const & rotset = rotamer_sets.has_rotamer_set_for_residue( other_seqpos ) ?
					rotamer_sets.rotamer_set_for_residue( other_seqpos )
					: fake_rotset;

				for ( Size irot = 0; irot <= rotset->num_rotamers(); irot++ ) {
					conformation::Residue const & res = irot == 0 ? pose.residue(other_seqpos) : *rotset->rotamer( irot );

					// Doing this because res.nbr_atom_xyz() causes crazy cache misses
					if ( ! is_protein || res.type().aa() == chemical::aa_gly ) {
						nbr_dist2 = res.nbr_atom_xyz().distance_squared( wat.water_xyz );
						if ( nbr_dist2 > longest_cb_water_interaction2 ) continue;
					}
					if ( nbr_dist2 > numeric::square(water_atom_range + res.nbr_radius() ) ) continue;

					LKD_ResidueInfo const & other_lkd_info = lkd_infos[ other_seqpos ][irot];


					Real occlude = 0;
					Real interact = 0;
					for ( Size iatom = 1; iatom <= res.nheavyatoms(); iatom++ ) {
						Real d2 = res.xyz(iatom).distance_squared( wat.water_xyz );
						if ( d2 > water_atom_range2 ) continue;

						Real lk_dome_iso_frac;
						Real lk_dome_frac;
						lk_dome_->single_water_atom_fractions( wat.base_xyz, wat.water_xyz,
							res.xyz(iatom), res.atom_type_index(iatom),
							lk_dome_iso_frac, lk_dome_frac );

						Size atom_offset = other_lkd_info.water_offset_for_atom()[iatom];

						Real bridge_frac = 0;
						if ( lk_dome_bridge_weight_ != 0 || lk_dome_bridge_uncpl_weight_ != 0 ) {
							for ( Size other_iwat = 1; other_iwat <= other_lkd_info.n_attached_waters()[iatom]; other_iwat++ ) {
								bridge_frac += lk_dome_->single_water_water_fraction_1way( wat.base_xyz, wat.water_xyz,
									other_lkd_info.waters()[other_iwat + atom_offset] );
							}
						}

						Real ball_bridge2_frac = 0;
						if ( lk_ball_bridge2_weight_ != 0 || lk_ball_bridge_uncpl2_weight_ != 0 ) {
							for ( Size other_iwat = 1; other_iwat <= other_lkd_info.n_attached_waters()[iatom]; other_iwat++ ) {
								ball_bridge2_frac += lk_dome_->single_water_water_bridge2_fraction_1way( wat.water_xyz,
									other_lkd_info.waters()[other_iwat + atom_offset] );
							}
						}

						// We're missing countpair here...

						Real final_score = wat.sol_value * (
							- lk_dome_iso_frac * lk_dome_iso_weight_
							- lk_dome_frac * lk_dome_weight_
							+ bridge_frac * lk_dome_bridge_weight_
							+ ball_bridge2_frac * lk_ball_bridge2_weight_ )
							+ bridge_frac * lk_dome_bridge_uncpl_weight_
							+ ball_bridge2_frac * lk_ball_bridge_uncpl2_weight_;


						Real this_occlude = lk_dome_->lk_ball()->get_lk_fractional_contribution_for_single_water(
							res.xyz(iatom), res.atom_type_index(iatom), wat.water_xyz );


						// if ( irot == 0 ) {
						//  bool found = false;
						//  for ( Size i = 0; i < water_assignments[water_seqpos][0].size(); i++) {
						//   if ( water_assignments[water_seqpos][0][i] == water_id ) {
						//    found = true;
						//    break;
						//   }
						//  }
						//  if ( found ) {
						//   std::cout << "H " << water_seqpos << " " << other_seqpos << " " << wat.iatom << " " << iatom << " " << "?" << " " << final_score << std::endl;
						//  }
						// }


						occlude += this_occlude;
						interact += final_score;

					}

					if ( occlude != 0 || interact != 0 ) {
						any_interactions = true;

						bit_set( irot, &this_mask[0] );
						bit_set( water_id, &oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos][irot][water_seqpos][0] );

						real_score_data.push_back( occlude );
						real_score_data.push_back( interact );

						largest_single_interaction_mag = std::max( largest_single_interaction_mag, std::abs( interact ) );
					}
				}

				if ( any_interactions ) {
					num_seqpos_added++;

					pos_water_pair_to_score_.push_back( (uint32_t) other_seqpos );
					pos_water_pair_to_score_.push_back( (uint32_t) (starting_size >> 1) * 3 + 1 ); // each pair gets expanded to 3 bytes
					for ( uint32_t bitfield : this_mask ) {
						pos_water_pair_to_score_.push_back( bitfield );
					}

				}
			}
			water_to_pwp_index_list_.push_back( num_seqpos_added );
		}
	}

	// Finish up riww

	seqpos_to_riww_at_pos_index_list_.clear();
	riww_at_pos_index_list_.clear();
	rotamer_interacts_with_waters_at_pos_.clear();

	// Loop over rotamer sets
	for ( Size other_seqpos = 1; other_seqpos <= pose.size(); other_seqpos++ ) {

		seqpos_to_riww_at_pos_index_list_.push_back( riww_at_pos_index_list_.size() + 1);

		// Loop over rotamers at this position
		Size irots_at_other = water_assignments[other_seqpos].size();
		for ( Size irot = 0; irot < irots_at_other; irot++ ) {
			riww_at_pos_index_list_.push_back( rotamer_interacts_with_waters_at_pos_.size() + 1);
			Size seqpos_added = 0;

			// Loop over positions this rotamer might interact with
			for ( Size water_seqpos = 1; water_seqpos <= pose.size(); water_seqpos++ ) {

				// Check if this rotamer interacts with this position
				bool all_zero = true;
				for ( uint32_t mask : oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos][irot][water_seqpos] ) {
					all_zero &= mask == 0;
				}
				if ( all_zero ) continue;

				// Interaction detected. Add the water_seqpos and masks to riww
				seqpos_added++;
				rotamer_interacts_with_waters_at_pos_.push_back( (uint32_t) water_seqpos );
				for ( uint32_t mask : oseqpos_orotamer_wseqpos_water_interact_mask[other_seqpos][irot][water_seqpos] ) {
					rotamer_interacts_with_waters_at_pos_.push_back( mask );
				}
			}
			riww_at_pos_index_list_.push_back( seqpos_added );
		}
	}

	// can store values 0 - 127
	//
	//



	// Ok, now we normalize the score_data

	int16_t max_occ_stored_value = 127;
	int32_t max_interact_stored_value = 32767;


	// This number should be occlusion_scale_ / FULLY_OCCLUDED;
	Real occlusion_scaler = lk_dome_->occlusion_max() / (FULLY_OCCLUDED + 0.9);

	Real occlusion_scale = lk_dome_->occlusion_max();

	occlusion_min_ = lk_dome_->occlusion_min() / occlusion_scale * FULLY_OCCLUDED;
	occlusion_span_ = FULLY_OCCLUDED - occlusion_min_;

	// Multiply the stored values by this to get the original value
	partial_score_scaler_ = largest_single_interaction_mag / (max_interact_stored_value + 0.9);
	// How to get the actual score
	score_scaler_ = partial_score_scaler_ / occlusion_span_;

	TR.Debug << "Max interaction magnitude: " << largest_single_interaction_mag << ". Score resolution is therefore: "
		<< partial_score_scaler_ << std::endl;


	Real inv_score_scaler = 1 / partial_score_scaler_;
	Real inv_occ_scaler = 1 / occlusion_scaler;


	score_data_.clear();
	score_data_.resize( (real_score_data.size() >> 1) * 3 );

	Size all_zero = 0;
	for ( Size i = 0; i < real_score_data.size() >> 1; i++ ) {

		score_data_[i*3+1] = int8_t(std::min<Real>( real_score_data[i*2+1], occlusion_scale ) * inv_occ_scaler);

		// Storing a 16 bit value into a 8 bit array.
		// If we always rely on the compiler to derefernce it, we don't have to worry about endianness
		int16_t store_interact = int16_t(real_score_data[i*2+2] * inv_score_scaler);
		*((int16_t *)&score_data_[i*3+2]) = store_interact;

		if ( score_data_[i*3+1] == 0 && score_data_[i*3+2] == 0 && score_data_[i*3+3] == 0 ) all_zero++;

		debug_assert( int64_t(std::min<Real>( real_score_data[i*2+1], occlusion_scale ) * inv_occ_scaler) <= max_occ_stored_value );
		debug_assert( int64_t(real_score_data[i*2+2] * inv_score_scaler) <= max_interact_stored_value );
	}

	TR.Debug << "Memory use: " << memory_use() << ". All zero: " << all_zero << " / " << (real_score_data.size()>>1) << std::endl;

}


// Set up all the stuff we need to perform the calculations from a fresh pose
core::pack::rotamer_set::RotamerSets
LK_DomeHelper::init_with_pose(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets_in
) {

	utility::vector1<utility::vector0<LKD_ResidueInfo>> lkd_infos;
	utility::vector1<utility::vector0<utility::vector0<Size>>> water_assignments;

	utility::pointer::shared_ptr< utility::vector1<utility::vector0<MyWaterHolder>> > waters =
		utility::pointer::make_shared< utility::vector1<utility::vector0<MyWaterHolder>> >();

	prepare_lkd_infos_and_assign_waters( pose, rotamer_sets_in, lkd_infos, *waters, water_assignments );

	prepare_simple_rotamer_values( *waters, water_assignments );

	prepare_score_data( pose, rotamer_sets_in, lkd_infos, *waters, water_assignments );

	resize_arrays( pose );

	if ( debug_store_waters_ ) {
		debug_waters_ = waters;
	}

	return rotamer_sets_in;
}


// This is an underestimate by about 1KB because we're only summing the big guys
Size
LK_DomeHelper::memory_use() const {

	Size size = 0;

	Size pose_size = seqpos_num_unique_waters_.size();
	Size total_waters = pose_size * waters_per_seqpos_;

	// The runtime stuff might not be populated yet

	size += pose_size * sizeof( decltype(current_rotamer_at_seqpos_)::value_type );
	size += pose_size * sizeof( decltype(current_score_at_seqpos_)::value_type );
	size += pose_size * sizeof( decltype(shadow_current_rotamer_at_seqpos_)::value_type );
	size += pose_size * sizeof( decltype(shadow_current_score_at_seqpos_)::value_type );
	size += pose_size * sizeof( decltype(shadow_mismatch_)::value_type );
	size += total_waters * sizeof( decltype(current_waters_)::value_type );
	size += total_waters * sizeof( decltype(current_water_occl_sum_)::value_type );
	size += total_waters * sizeof( decltype(current_water_interact_sum_)::value_type );
	size += total_waters * sizeof( decltype(shadow_current_waters_)::value_type );
	size += total_waters * sizeof( decltype(shadow_current_water_occl_sum_)::value_type );
	size += total_waters * sizeof( decltype(shadow_current_water_interact_sum_)::value_type );
	size += seqpos_to_rotamer_to_waters_masks_.size() * sizeof( decltype(seqpos_to_rotamer_to_waters_masks_)::value_type );
	size += seqpos_to_water_to_pwp_index_list_.size() * sizeof( decltype(seqpos_to_water_to_pwp_index_list_)::value_type );
	size += seqpos_to_riww_at_pos_index_list_.size() * sizeof( decltype(seqpos_to_riww_at_pos_index_list_)::value_type );
	size += seqpos_num_unique_waters_.size() * sizeof( decltype(seqpos_num_unique_waters_)::value_type );
	size += seqpos_num_rotamers_32_.size() * sizeof( decltype(seqpos_num_rotamers_32_)::value_type );
	size += rotamer_to_waters_masks_.size() * sizeof( decltype(rotamer_to_waters_masks_)::value_type );
	size += water_to_pwp_index_list_.size() * sizeof( decltype(water_to_pwp_index_list_)::value_type );
	size += pos_water_pair_to_score_.size() * sizeof( decltype(pos_water_pair_to_score_)::value_type );
	size += score_data_.size() * sizeof( decltype(score_data_)::value_type );
	size += riww_at_pos_index_list_.size() * sizeof( decltype(riww_at_pos_index_list_)::value_type );
	size += rotamer_interacts_with_waters_at_pos_.size() * sizeof( decltype(rotamer_interacts_with_waters_at_pos_)::value_type );

	return size;

}


core::scoring::lkball::LK_DomeEnergyCOP
LK_DomeHelper::lk_dome() const { return lk_dome_; }

Real
LK_DomeHelper::current_score() const {
	return shadow_total_score_ * score_scaler_;
}


}
} //guidance_scoreterms
} //pack
} //core
