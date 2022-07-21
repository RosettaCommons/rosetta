// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/guidance_scoreterms/sap/SapConstraintHelper.cc
/// @brief  Score term that applies the RotamerPSSMConstraint as an energy
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>
#include <core/pack/guidance_scoreterms/sap/util.hh>


// Package headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/init_id_map.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/conformation/util.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/select/util.hh>

// Basic Headers
#include <basic/Tracer.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <utility/graph/Graph.hh>

#include <boost/format.hpp>

#include <core/conformation/symmetry/SymmetryInfo.hh> // AUTO IWYU For SymmetryInfo
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask
#include <core/scoring/methods/EnergyMethodOptions.hh> // AUTO IWYU For EnergyMethodOptions
#include <core/select/residue_selector/ResidueSelector.hh> // AUTO IWYU For ResidueSelector
#include <utility/stream_util.hh> // AUTO IWYU For operator<<


static basic::Tracer TR("core.pack.guidance_scoreterms.sap.SapConstraintHelper");

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {




// Here is a brief description of the three operating modes:

// Slow: This performs nearly the full sap calculation except we use an approximation to sasa.
//         Sasa is calculated by the "blocks" system where fractions of an atom's sasa are "blocked"
//          by neighboring atoms. There are then sigmoids fit to each atom type that give a sasa value
//         The sasa's are used to calculate SASA scores which are what are actually stored in the sasa_scores
//          these are the values that can be summed in each 5A window to produce the sap scores
//         The sap values are the final array which just keep track of the sum in each atom's 5A window
//         Summing the sap values that are > 0 gives the sap score

// Fast: This is just like slow except the Sasa values are pre-computed for each rotamer. During the pre-computation
//         we assume that any residue that can design only has a backbone, CA, and CB atom

// Lightning: The biggest change going into lightning is that we assume that all rotamers of a given AA-type behave the same
//         Additionally, we don't do a per-atom 0-clip anymore but instead sum up all the atoms on a rotamer and then clip to 0




SapConstraintHelper::SapConstraintHelper( SapConstraintOptionsCOP const & options ) :
	options_( options )
{
	init();
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// PACKER RUNTIME ROUTINES /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// These functions have to be as fast as possible as they determine the packer slowdown
// High mem-use is almost always preferable to code speed
//   -- But don't forget about the cache! Sometimes smaller is better in that regard


// This is the main function that gets called during packing. Notably this function either
//  adjusts the current internal resvector or restarts from scratch.
//  Before every substition, the shadow is copied back into our working registers so that
//  we can easily recover from failed packer moves.
core::Real
SapConstraintHelper::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	core::Size const substitution_position
) {

	if ( symm_info_ ) return symm_calculate_energy( resvect, substitution_position );

	if ( substitution_position == 0 ) {
		reinit_with_resvect( resvect );
	} else {
		restore_from_shadow();

		// If it's not a sasa position, it's never important
		if ( ! sasa_positions_[ substitution_position ] ) return current_score_;

		// If it's not a calculate_position, it only matters in "slow" where sasa are changing
		if ( fast_ && ! sap_calculate_positions_[ substitution_position ] ) return current_score_;

		recalc_positions_scratch_[1] = substitution_position;

		if ( lightning_ ) {
			if ( internal_resvect_[substitution_position]->aa() != resvect[substitution_position]->aa() ) {
				add_remove_rotamer_lightning( internal_resvect_, substitution_position, false );
				add_remove_rotamer_lightning( resvect, substitution_position, true );
			}
		} else {
			if ( fast_ ) {
				add_remove_rotamer_fast( internal_resvect_, substitution_position, false );
				add_remove_rotamer_fast( resvect, substitution_position, true );
			} else {
				add_remove_rotamer( internal_resvect_, substitution_position, false );
				add_remove_rotamer( resvect, substitution_position, true );
				recalculate_sasa();
			}

			recalculate_saps( recalc_positions_scratch_ );
		}
	}

	return current_score_;
}

// All we need to do to commit a substitution is to save our current registers into the shadow
void
SapConstraintHelper::commit_considered_substitution() {
	// std::cout << "   commit" << std::endl;
	save_to_shadow();
}

// This function copies all of the shadow registers back into the working registers. This is used
//  so that we can ignore a failed packer substitution
void
SapConstraintHelper::restore_from_shadow() {

	for ( Size seqpos = 1; seqpos <= internal_resvect_.size(); seqpos ++ ) {
		if ( ! shadow_mismatch_[seqpos] ) continue;
		shadow_mismatch_[seqpos] = false;

		if ( lightning_ ) {
			lightning_current_res_sap_[seqpos] = lightning_current_res_sap_shadow_[seqpos];
			internal_resvect_[seqpos] = internal_resvect_shadow_[seqpos];
		} else {
			atom_sasa_score_[seqpos] = atom_sasa_score_shadow_[seqpos];
			sasa_blocks_[seqpos] = sasa_blocks_shadow_[seqpos];
			atom_sap_[seqpos] = atom_sap_shadow_[seqpos];
			block_param_offset_[seqpos] = block_param_offset_shadow_[seqpos];
			internal_resvect_[seqpos] = internal_resvect_shadow_[seqpos];
			atom_sasa_score_fast_[seqpos] = atom_sasa_score_fast_shadow_[seqpos];
		}
		if ( symm_info_ ) {
			symm_work_resvect_[seqpos] = symm_work_resvect_shadow_[seqpos];
		}

	}

	current_score_ = current_score_shadow_;
}

// Save our corrent working registers into the shadow. This is what it means to "commit" a change.
void
SapConstraintHelper::save_to_shadow() {

	for ( Size seqpos = 1; seqpos <= internal_resvect_.size(); seqpos ++ ) {
		if ( ! shadow_mismatch_[seqpos] ) continue;
		shadow_mismatch_[seqpos] = false;

		if ( lightning_ ) {
			lightning_current_res_sap_shadow_[seqpos] = lightning_current_res_sap_[seqpos];
			internal_resvect_shadow_[seqpos] = internal_resvect_[seqpos];
		} else {
			atom_sasa_score_shadow_[seqpos] = atom_sasa_score_[seqpos];
			sasa_blocks_shadow_[seqpos] = sasa_blocks_[seqpos];
			atom_sap_shadow_[seqpos] = atom_sap_[seqpos];
			block_param_offset_shadow_[seqpos] = block_param_offset_[seqpos];
			internal_resvect_shadow_[seqpos] = internal_resvect_[seqpos];
			atom_sasa_score_fast_shadow_[seqpos] = atom_sasa_score_fast_[seqpos];
		}
		if ( symm_info_ ) {
			symm_work_resvect_shadow_[seqpos] = symm_work_resvect_[seqpos];
		}

	}

	current_score_shadow_ = current_score_;
}


// We're handling symmetry in a non-symmetric way
//  This function is basically an exact copy of calculate_energy(), except we go through the remove/add
//   process for each asu individually. Then, since the score_positions_ only refers to
//   the asu, we return only the asu score
core::Real
SapConstraintHelper::symm_calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	core::Size const substitution_position
) {

	if ( substitution_position == 0 ) {
		reset_calculation(); // Need this here to resize symm_work_resvect_

		// All we need to do is fill symm_work_resvect_ with our rotamers and then call calculate_energy again
		for ( Size seqpos = 1; seqpos <= resvect.size(); seqpos++ ) {
			if ( symm_info_->bb_follows( seqpos ) != 0 ) continue; // not asu

			// store the asu rotamer because we actually used those
			symm_work_resvect_[ seqpos ] = resvect[ seqpos ];

			// store the other rotamers that we generated earlier
			utility::vector1< core::conformation::ResidueCOP > const & other_rots = symm_rotamer_to_other_rotamers_.at( &*resvect[ seqpos ] );
			for ( core::conformation::ResidueCOP const & other_rot : other_rots ) {
				symm_work_resvect_[ other_rot->seqpos() ] = other_rot;
			}
		}

		// This has a vector of all known residues, so this function just works
		reinit_with_resvect( symm_work_resvect_, true );
	} else {
		debug_assert( symm_info_->bb_follows( substitution_position ) == 0 );

		restore_from_shadow();

		// If it's not a sasa position, it's never important
		if ( ! sasa_positions_[ substitution_position ] ) return current_score_;

		// If it's not a calculate_position, it only matters in "slow" where sasa are changing
		if ( fast_ && ! sap_calculate_positions_[ substitution_position ] ) return current_score_;

		// Get all the symmetric copy rotamers
		utility::vector1< core::conformation::ResidueCOP > const & other_rots =
			symm_rotamer_to_other_rotamers_.at( &*resvect[ substitution_position ] );

		utility::vector1<Size> recalc_positions;

		// Loop over asu rotamer and symmetric copy rotamers
		for ( Size irot = 0; irot <= other_rots.size(); irot++ ) {
			core::conformation::ResidueCOP const & rotamer = irot == 0 ? resvect[ substitution_position ] : other_rots[ irot ];
			Size seqpos = rotamer->seqpos();

			// add one by one to the symm_work_resvect_ and mark a mismatch
			symm_work_resvect_[ seqpos ] = rotamer;
			shadow_mismatch_[ seqpos ] = true;

			// Exactly the same as calculate_energy() except we defer the recalculate_ steps
			if ( lightning_ ) {
				if ( internal_resvect_[seqpos]->aa() != resvect[seqpos]->aa() ) {
					add_remove_rotamer_lightning( internal_resvect_, seqpos, false );
					add_remove_rotamer_lightning( symm_work_resvect_, seqpos, true );
				}
			} else {
				if ( fast_ ) {
					add_remove_rotamer_fast( internal_resvect_, seqpos, false );
					add_remove_rotamer_fast( symm_work_resvect_, seqpos, true );
				} else {
					add_remove_rotamer( internal_resvect_, seqpos, false );
					add_remove_rotamer( symm_work_resvect_, seqpos, true );
				}
				recalc_positions.push_back( seqpos );
			}
		}

		if ( ! fast_ ) recalculate_sasa();
		if ( ! lightning_ ) recalculate_saps( recalc_positions );

	}
	return current_score_;
}




// Recalculate the sap scores of all positions in positions_to_update. Notably, the atom_sasa_scores
//  need to be correct before we can call this.
//  As this function loops through each sequence position, it first clears out the sap score from
//  all the atoms that are currently there.
//  Then, for each position that interacts with this one it looks to see if there are known interactions
//  between the rotamer currently there and us. (Interaction here is atoms within 5Å)
//  Finally, if there is an interaction, we loop through all interacting pairs of atoms and accumulate
//  the atom scores from the other residue to us
//  This function is only used in "slow" and "fast"
void
SapConstraintHelper::recalculate_saps( utility::vector1<Size> const & positions_to_update ) {

	for ( Size seqpos : positions_to_update ) {
		if ( !score_positions_[seqpos] ) continue;

		core::conformation::ResidueCOP const & internal_res = internal_resvect_[seqpos];
		Size first_sidechain = internal_res->first_sidechain_atom();
		Size natoms = my_natoms( internal_res );

		// Clear out current sap scores for this position
		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain ; iat++ ) {
			float sap = atom_sap_[seqpos][iat];
			if ( sap > 0 ) current_score_ -= sap;
			atom_sap_[seqpos][iat] = 0;
		}

		utility::vector1< Size > check_positions = check_positions_sap_[ seqpos ];

		std::pair< conformation::Residue const *, conformation::Residue const * > key(
			&*internal_resvect_[seqpos], &*internal_resvect_[seqpos] );

		// Loop through all other positions known to interact with this one
		for ( Size position : check_positions_sap_[seqpos] ) {

			key.second = &*internal_resvect_[position];
			auto iter = atom_within_5_map_.find( key );

			// The current rotamer pair doesn't have any interactions
			if ( iter == atom_within_5_map_.end() ) continue;

			// In fast_, the atom_sasa_scores don't change so we just grab the pointer to the precomputed
			float * other_sasa_scores = fast_ ? atom_sasa_score_fast_[position] : &atom_sasa_score_[position].front();

			// atom_within_5 tells you which atoms on key.second are within 5Å of each atom on key.first
			// Remember: This is the structure of atom_within_5_value
			//           std::pair< Size, uint8_t[ATOM_WITHIN_5_ELEMS]>
			// The uint8 has two sections:
			//   1. A list of length ( natoms - first_sidechain + 1 ) that contains the number of key.second atoms for each atom
			//   2. A bunch of lists of key.second atoms (one list for each atom on key.first but all concated together)
			//
			// There are only 40 bytes to store 2. however. Once it overflows, the Size refers to the offset into
			//   atom_within_5_ where the list continues

			// Get atom_within_5_value
			atom_within_5_value const & value = iter->second;

			// offset is offset into atom_within_5_
			// local offset is offset into value.second
			Size offset = value.first;
			Size local_offset = natoms - first_sidechain + 1;

			Size iat = 0;   // which iat on key.first
			Size elems = 0; // how many atoms on key.second for this iat
			Size ielem = 0; // which atom in elems are we on

			// Here we iterate inside value.second;
			for ( ; (int)iat <= (int)natoms - (int)first_sidechain && local_offset < ATOM_WITHIN_5_ELEMS; iat++ ) {

				elems = value.second[iat];
				for ( ielem = 0; ielem < elems && local_offset < ATOM_WITHIN_5_ELEMS; ielem++, local_offset++ ) {

					// Accumulate sasa_score from other atom onto our atom
					Size other_atom = value.second[local_offset];
					atom_sap_[seqpos][iat] += other_sasa_scores[other_atom];
				}
			}

			// This is a no-op here (it's always +0 )but keeping it for consistency with other sections
			offset += local_offset - ATOM_WITHIN_5_ELEMS;

			// Finish up the current atom if the jump from atom_within_5_value to atom_within_5_ occurred mid-atom
			for ( ; ielem < elems; ielem++, offset++ ) {

				// Accumulate sasa_score from other atom onto our atom
				Size other_atom = atom_within_5_[offset];
				atom_sap_[seqpos][iat-1] += other_sasa_scores[other_atom];
			}

			// Here we iterate inside atom_within_5_
			for ( ; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {

				elems = value.second[iat];
				for ( ielem = 0; ielem < elems; ielem++, offset++ ) {

					// Accumulate sasa_score from other atom onto our atom
					Size other_atom = atom_within_5_[offset];
					atom_sap_[seqpos][iat] += other_sasa_scores[other_atom];
				}
			}
		}


		// Now we accumulate our atom_sap_ into current_score_
		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain ; iat++ ) {
			float sap = atom_sap_[seqpos][iat];
			if ( sap > 0 ) current_score_ += sap;
		}
	}
}


// This function recalculates the sasa_scores from the sasa_blocks_.
//  Since we approximate sasa by using sasa_blocks_, this function updates
//  any sasa_scores that have gone stale.
//  This function loops through all sequence positions and then for any that
//  have dirty atoms, it simply uses the BlockParams to calculate the current
//  sasa score from the current sasa_blocks_
//  After updating the sasa_scores for the current rotamer, it calls update_neighbors_sap()
//  so that any other residues that depend on this one get their sap scores updated
//  This function is only used in "slow"

// This function clears the dirty buffer
void
SapConstraintHelper::recalculate_sasa() {

	// Only residues with shadow_mismatch_ will have any dirty atoms
	for ( Size seqpos = 1; seqpos <= internal_resvect_.size(); seqpos++ ) {
		if ( ! shadow_mismatch_[seqpos] ) continue;

		core::conformation::ResidueCOP const & internal_res = internal_resvect_[seqpos];
		Size first_sidechain = internal_res->first_sidechain_atom();
		Size natoms = my_natoms( internal_res );

		Size block_param_offset = block_param_offset_[ seqpos ];

		bool any_delta_sasa = false;

		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain ; iat++ ) {
			if ( ! dirty_sasa_[seqpos][iat] ) continue;

			BlockParam const & param = all_block_params_[block_param_offset+iat];

			float sasa;
			uint16_t blocks = sasa_blocks_[seqpos][iat];
			float frac_sasa = 0;
			if ( blocks > param.full_block ) {
				sasa = 0;
			} else if ( blocks < param.no_block ) {
				sasa = param.max_sasa_score;
				//frac_sasa = 1;
			} else {
				frac_sasa = float( param.full_block - blocks ) / float( param.full_block - param.no_block );
				sasa = frac_sasa * param.max_sasa_score;
			}

			float delta = sasa - atom_sasa_score_[seqpos][iat];
			delta_sasa_scratch_[iat] = delta;
			any_delta_sasa |= delta != 0;

			atom_sasa_score_[seqpos][iat] = sasa;
		}

		if ( any_delta_sasa ) {
			update_neighbors_sap( &delta_sasa_scratch_.front(), seqpos, false, false );
		}

		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain ; iat++ ) {
			dirty_sasa_[seqpos][iat] = false;
		}
	}
}

// This function updates all sap neighbors of a residue when its sasa values change.
//  The first thing to do is to loop through all other positions that are known to have
//  rotamers that interact with this position. Then, for each of those positions, we loop
//  through all atoms pairs using the atom_within_5 and make the necessary adjustments to
//  the sap values
//  This function is only used by "slow" and "lightning"

// Called when a rotamer changes its sasa values and needs to update all the saps relying on it
// Only checks dirty atom positions
void
SapConstraintHelper::update_neighbors_sap(
	float const * delta_sasa_score,
	Size seqpos,
	bool invert,
	bool skip_self
) {
	float multiplier = invert ? -1 : 1;

	core::conformation::ResidueCOP const & internal_res = internal_resvect_[seqpos];
	Size first_sidechain = internal_res->first_sidechain_atom();
	Size natoms = my_natoms( internal_res );

	utility::vector1< Size > check_positions = check_positions_sap_[ seqpos ];

	std::pair< conformation::Residue const *, conformation::Residue const * > key(
		&*internal_resvect_[seqpos], &*internal_resvect_[seqpos] );

	for ( Size position : check_positions_sap_[seqpos] ) {
		if ( skip_self && position == seqpos ) continue;
		if ( ! score_positions_[position] ) continue;

		key.second = &*internal_resvect_[position];
		auto iter = atom_within_5_map_.find( key );

		if ( iter == atom_within_5_map_.end() ) continue;

		bool created_shadow_mismatch = false;

		// atom_within_5 tells you which atoms on key.second are within 5Å of each atom on key.first
		// Remember: This is the structure of atom_within_5_value
		//           std::pair< Size, uint8_t[ATOM_WITHIN_5_ELEMS]>
		// The uint8 has two sections:
		//   1. A list of length ( natoms - first_sidechain + 1 ) that contains the number of key.second atoms for each atom
		//   2. A bunch of lists of key.second atoms (one list for each atom on key.first but all concated together)
		//
		// There are only 40 bytes to store 2. however. Once it overflows, the Size refers to the offset into
		//   atom_within_5_ where the list continues

		// Get atom_within_5_value
		atom_within_5_value const & value = iter->second;

		// offset is offset into atom_within_5_
		// local offset is offset into value.second
		Size offset = value.first;
		Size local_offset = natoms - first_sidechain + 1;

		Size iat = 0;   // which iat on key.first
		Size elems = 0; // how many atoms on key.second for this iat
		Size ielem = 0; // which atom in elems are we on
		float delta_amount = 0;

		// Here we iterate inside value.second;
		for ( ; (int)iat <= (int)natoms - (int)first_sidechain && local_offset < ATOM_WITHIN_5_ELEMS; iat++ ) {

			elems = value.second[iat];
			if ( elems > 0 ) {
				if ( dirty_sasa_[seqpos][iat] ) {
					created_shadow_mismatch = true;

					delta_amount = delta_sasa_score[iat] * multiplier;

					for ( ielem = 0; ielem < elems && local_offset < ATOM_WITHIN_5_ELEMS; ielem++, local_offset++ ) {

						Size other_atom = value.second[local_offset];

						// Accumulate the delta into the other position
						add_to_sap( delta_amount, position, other_atom );
					}
				} else {
					local_offset += elems;
					elems = 0;
				}
			}
		}

		// If we went past the end of the local buffer (because of the local_offset += elems line), the correct for that
		offset += local_offset - ATOM_WITHIN_5_ELEMS;

		// Finish up the current atom if the jump from atom_within_5_value to atom_within_5_ occurred mid-atom
		for ( ; ielem < elems; ielem++, offset++ ) {

			Size other_atom = atom_within_5_[offset];

			// Accumulate the delta into the other position
			add_to_sap( delta_amount, position, other_atom );
		}

		// Here we iterate inside atom_within_5_
		for ( ; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {

			elems = value.second[iat];
			if ( elems > 0 ) {
				if ( dirty_sasa_[seqpos][iat] ) {
					created_shadow_mismatch = true;

					delta_amount = delta_sasa_score[iat] * multiplier;

					for ( ielem = 0; ielem < elems; ielem++, offset++ ) {

						Size other_atom = atom_within_5_[offset];

						// Accumulate the delta into the other position
						add_to_sap( delta_amount, position, other_atom );
					}
				} else {
					offset+= elems;
				}
			}
		}

		// If we actually changed something, the shadow is out of wack
		if ( created_shadow_mismatch ) {
			shadow_mismatch_[position] = true;
		}
	}
}

// Don't add any fancy logic here
// Purely a helper function to deal with +- around 0
void
SapConstraintHelper::add_to_sap( float amount, Size seqpos, Size iat ) {
	debug_assert( score_positions_[seqpos] );   // if this gets called we're wasting cpu time

	debug_assert( iat + internal_resvect_[seqpos]->first_sidechain_atom() <= internal_resvect_[seqpos]->natoms()
		&& ! internal_resvect_[seqpos]->atom_is_backbone( iat + internal_resvect_[seqpos]->first_sidechain_atom() ) );

	if ( amount == 0 ) return;


	float orig_sap = atom_sap_[seqpos][iat];
	float new_sap = orig_sap + amount;
	atom_sap_[seqpos][iat] = new_sap;

	if ( amount > 0 ) {
		if ( orig_sap >= 0 ) {
			current_score_ += amount;   // was positive the whole time
		} else {
			if ( new_sap > 0 ) {
				current_score_ += new_sap; // we made it positive
			}
		}
	} else {
		if ( orig_sap > 0 ) {
			if ( new_sap < 0 ) {
				current_score_ -= orig_sap; // me made it negative
			} else {
				current_score_ += amount; // was positive the whole time
			}
		}
	}
}

// Don't add any fancy logic here
// Purely a helper function to deal with +- around 0
void
SapConstraintHelper::lightning_add_to_sap( float amount, Size seqpos  ) {
	debug_assert( score_positions_[seqpos] );   // if this gets called we're wasting cpu time

	if ( amount == 0 ) return;

	float orig_sap = lightning_current_res_sap_[seqpos];
	float new_sap = orig_sap + amount;
	lightning_current_res_sap_[seqpos] = new_sap;

	if ( amount > 0 ) {
		if ( orig_sap >= 0 ) {
			current_score_ += amount;   // was positive the whole time
		} else {
			if ( new_sap > 0 ) {
				current_score_ += new_sap; // we made it positive
			}
		}
	} else {
		if ( orig_sap > 0 ) {
			if ( new_sap < 0 ) {
				current_score_ -= orig_sap; // me made it negative
			} else {
				current_score_ += amount; // was positive the whole time
			}
		}
	}
}

// Accumulate the sasa_scores from the current rotamer at seqpos into all neighboring
//  positions. Since in lightning there are no per-atom scores, this function only needs to lookup
//  a single value for each residue-residue pair

// skip self by default
void
SapConstraintHelper::lightning_update_neighbors_sap(
	Size seqpos,
	bool invert,
	bool dont_update_self
) {
	float multiplier = invert ? -1 : 1;

	core::conformation::ResidueCOP const & internal_res = internal_resvect_[seqpos];

	// std::cout << internal_res->name() << std::endl;

	utility::vector1< Size > check_positions = check_positions_sap_[ seqpos ];

	for ( Size position : check_positions_sap_[seqpos] ) {
		if ( position == seqpos ) continue;
		if ( ! score_positions_[position] ) continue;


		core::conformation::ResidueCOP const & other_res = internal_resvect_[position];
		if ( &*other_res == &*fake_rotamer_ ) continue;

		if ( seqpos < position ) {

			std::pair< float, float > const & on_seq_on_pos = lightning_2b_lookup( internal_res->aa(), seqpos,
				other_res->aa(), position );

			if ( on_seq_on_pos.first == 0 && on_seq_on_pos.second  == 0 ) continue;

			shadow_mismatch_[ position ] = true;

			if ( ! dont_update_self ) {
				lightning_add_to_sap( on_seq_on_pos.first * multiplier, seqpos );
			}
			lightning_add_to_sap( on_seq_on_pos.second * multiplier, position );
		} else {
			std::pair< float, float > const & on_pos_on_seq = lightning_2b_lookup( other_res->aa(), position,
				internal_res->aa(), seqpos );

			if ( on_pos_on_seq.first == 0 && on_pos_on_seq.second  == 0 ) continue;

			shadow_mismatch_[ position ] = true;

			lightning_add_to_sap( on_pos_on_seq.first * multiplier, position );
			if ( ! dont_update_self ) {
				lightning_add_to_sap( on_pos_on_seq.second * multiplier, seqpos );
			}
		}
	}
}

// Add or remove a rotamer during lightning packing
// If we're adding, we need to put the new rotamer into the internal_resvect
//  and then assign it a starting sap of it's 1b sap energy
// If we are removing, we simply clear out the current sap score
// Finally, we update all the neighbors to account for the change (and if we
//  are adding, that function also accumulates the 2b energies into our new residue)

// If removing rotamers, this function always expects
// to immediately be called again with adding a rotamer to the same positions
void
SapConstraintHelper::add_remove_rotamer_lightning(
	utility::vector1< core::conformation::ResidueCOP > const & resvect,
	Size const substitution_position,
	bool add
) {
	shadow_mismatch_[substitution_position] = true;

	if ( add ) {
		core::conformation::ResidueCOP const & rotamer = resvect[ substitution_position ];
		internal_resvect_[substitution_position] = resvect[ substitution_position ];

		float self_sap = lightning_1b_lookup( rotamer->aa(), substitution_position );
		if ( self_sap > 0 ) current_score_ += self_sap;

		debug_assert( lightning_current_res_sap_[ substitution_position ] == 0 );
		lightning_current_res_sap_[ substitution_position ] = self_sap;

	} else {
		if ( score_positions_[substitution_position] ) {

			float cur_sap = lightning_current_res_sap_[ substitution_position ];
			if ( cur_sap > 0 ) {
				current_score_ -= cur_sap;
			}
			lightning_current_res_sap_[ substitution_position ] = 0;
		}
	}

	lightning_update_neighbors_sap( substitution_position, ! add, ! add  );
}

// Add or remove a rotamer during fast packing
// If we're adding, we need to put the new rotamer into the internal_resvect
//  and then swap out the current sasa vector
// If we are removing, we simply clear out the current sap score
// Finally, we update all the neighbors to account for the change (and if we
//  are adding, that function also accumulates sap into our new residue)

// If removing rotamers, this function always expects
// to immediately be called again with adding a rotamer to the same positions
void
SapConstraintHelper::add_remove_rotamer_fast(
	utility::vector1< core::conformation::ResidueCOP > const & resvect,
	Size const substitution_position,
	bool add
) {
	shadow_mismatch_[substitution_position] = true;

	if ( add ) {
		// Reset sasa_blocks_ and dirty_sasa_ for this rotamer
		core::conformation::ResidueCOP const & rotamer = resvect[ substitution_position ];

		internal_resvect_[substitution_position] = resvect[ substitution_position ];
		atom_sasa_score_fast_[substitution_position] = rotamer_to_sasa_data_.at( &*rotamer );

	} else {

		if ( score_positions_[substitution_position] ) {
			core::conformation::ResidueCOP const & prev_rotamer = resvect[ substitution_position ];
			Size first_sidechain = prev_rotamer->first_sidechain_atom();
			Size natoms = my_natoms( prev_rotamer );

			for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {

				float sap = atom_sap_[ substitution_position ][ iat ];
				if ( sap > 0 ) current_score_ -= sap;
				atom_sap_[ substitution_position ][ iat ] = 0;
			}
		}
	}


	update_neighbors_sap( atom_sasa_score_fast_[ substitution_position ], substitution_position, ! add, true );

}
// Add or remove a rotamer during slow packing
// If we're adding, we need to put the new rotamer into the internal_resvect
//  and clear out the sasa_blocks and set all its atoms to dirty
// If we are removing, we first removing the sap contribution from all neighboring
//  residues caused by this residue and then clear out the saps and sasa_scores
// Finally, in both cases, we need to go through and update all of the sasa_block values
//  for all of the neighboring residues

// If removing rotamers, this function always expects
// to immediately be called again with adding a rotamer to the same positions
void
SapConstraintHelper::add_remove_rotamer(
	utility::vector1< core::conformation::ResidueCOP > const & resvect,
	Size const substitution_position,
	bool add
) {

	shadow_mismatch_[substitution_position] = true;

	if ( add ) {
		// Reset sasa_blocks_ and dirty_sasa_ for this rotamer
		core::conformation::ResidueCOP const & rotamer = resvect[ substitution_position ];
		Size first_sidechain = rotamer->first_sidechain_atom();
		Size natoms = my_natoms( rotamer );

		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
			dirty_sasa_[substitution_position][iat] = true;
			sasa_blocks_[ substitution_position ][ iat ] = 0;

		}

		internal_resvect_[substitution_position] = resvect[ substitution_position ];
		block_param_offset_[substitution_position] = rotamer_to_block_param_offset_.at( &*rotamer );

	} else {

		// Remove score info for previous rotamer

		core::conformation::ResidueCOP const & prev_rotamer = resvect[ substitution_position ];
		Size first_sidechain = prev_rotamer->first_sidechain_atom();
		Size natoms = my_natoms( prev_rotamer );

		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
			dirty_sasa_[substitution_position][iat] = true;
		}

		// First clear out all saps that rely on these sasa
		update_neighbors_sap( &atom_sasa_score_[ substitution_position ].front(), substitution_position, true, true );

		// Next clear out all of this rotamers saps and sasas
		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {

			if ( score_positions_[substitution_position] ) {
				float sap = atom_sap_[ substitution_position ][ iat ];
				if ( sap > 0 ) current_score_ -= sap;
				atom_sap_[ substitution_position ][ iat ] = 0;
			}

			atom_sasa_score_[ substitution_position ][ iat ] = 0;
			dirty_sasa_[substitution_position][iat] = false;
		}

	}

	// All positions that might have sasa interactions from us
	utility::vector1< Size > check_positions = check_positions_block_[ substitution_position ];

	std::pair< conformation::Residue const *, conformation::Residue const * > keyl(
		&*resvect[substitution_position], &*resvect[substitution_position] );
	std::pair< conformation::Residue const *, conformation::Residue const * > keyh(
		&*resvect[substitution_position], &*resvect[substitution_position] );

	// loop through all potentially interacting positions and for positions
	// where there are interactions, add to blocks and mark sap dirty

	for ( Size position : check_positions_block_[substitution_position] ) {

		bool is_before_us = position < substitution_position;


		Size offset = 0;
		if ( is_before_us ) {
			keyl.first = &*resvect[position];
			auto iter = interacting_block_offset_.find( keyl );
			if ( iter == interacting_block_offset_.end() ) continue;
			offset = iter->second;
		} else {
			keyh.second = &*resvect[position];
			auto iter = interacting_block_offset_.find( keyh );
			if ( iter == interacting_block_offset_.end() ) continue;
			offset = iter->second;
		}


		shadow_mismatch_[position] = true;


		if ( is_before_us ) {
			store_sasa_blocks( offset, position, add, true );
			store_sasa_blocks( offset, substitution_position, add, false );

		} else {
			store_sasa_blocks( offset, substitution_position, add, false );
			store_sasa_blocks( offset, position, add, true );
		}
	}
}


// Each rotamer pair knows how many blocks it causes on each of the other residue's atoms
//  These have been accumulated ahead of time, so all we need to do is loop through each atom
//  and apply the pre-computed block value
//  This function is only called in "slow"
void
SapConstraintHelper::store_sasa_blocks( Size & offset, Size seqpos, bool add, bool mark_dirty ) {

	core::conformation::ResidueCOP const & rotamer = internal_resvect_[ seqpos ];
	Size first_sidechain = rotamer->first_sidechain_atom();
	Size natoms = my_natoms( rotamer );

	for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
		uint8_t block = interacting_block_[offset++];

		if ( add ) {
			sasa_blocks_[seqpos][iat] += block;
		} else {
			debug_assert( sasa_blocks_[seqpos][iat] >= block );
			sasa_blocks_[seqpos][iat] -= block;
		}

		if ( mark_dirty && block > 0 ) {
			dirty_sasa_[seqpos][iat] = true;
		}
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
SapConstraintHelper::reinit_with_resvect(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	bool skip_reset /* = false */
) {
	if ( ! skip_reset ) reset_calculation();

	// By using fake_rotamer_ sasas will only be calculated for stuff that already exists
	utility::vector1< core::conformation::ResidueCOP > in_progress( resvect.size(), fake_rotamer_ );
	utility::vector1<Size> recalc_positions;

	for ( Size seqpos = 1; seqpos <= resvect.size(); seqpos++ ) {
		// If it's not a sasa position, it's never important
		if ( ! sasa_positions_[ seqpos ] ) continue;

		// If it's not a calculate_position, it only matters in "slow" where sasa are changing
		if ( fast_ && ! sap_calculate_positions_[ seqpos ] ) continue;

		in_progress[seqpos] = resvect[seqpos];
		if ( lightning_ ) {
			add_remove_rotamer_lightning( in_progress, seqpos, true );
		} else if ( fast_ ) {
			add_remove_rotamer_fast( in_progress, seqpos, true );
		} else {
			add_remove_rotamer( in_progress, seqpos, true );
		}
		recalc_positions.push_back(seqpos);
	}

	if ( ! lightning_ ) {
		if ( ! fast_ ) {
			recalculate_sasa();
		}
		recalculate_saps( recalc_positions );
	}

	save_to_shadow();
}


// Clear and resize all arrays that accumulate things related to scores.
// This will not prepare for a new pose however
void
SapConstraintHelper::reset_calculation() {

	Size pose_size = score_positions_.size();

	internal_resvect_.clear();
	internal_resvect_.resize( pose_size, fake_rotamer_ );

	block_param_offset_.clear();
	block_param_offset_.resize( pose_size );
	sasa_blocks_.clear();
	sasa_blocks_.resize( pose_size, utility::vector0<uint8_t>(max_rotamer_atoms_) );
	atom_sasa_score_.clear();
	atom_sasa_score_.resize( pose_size, utility::vector0<float>(max_rotamer_atoms_) );
	atom_sasa_score_fast_.clear();
	atom_sasa_score_fast_.resize( pose_size, &delta_sasa_zeros_.front() );
	atom_sap_.clear();
	atom_sap_.resize( pose_size, utility::vector0<float>(max_rotamer_atoms_) );

	current_score_ = 0;

	internal_resvect_shadow_.clear();
	internal_resvect_shadow_.resize( pose_size, fake_rotamer_ );

	block_param_offset_shadow_.clear();
	block_param_offset_shadow_.resize( pose_size );
	sasa_blocks_shadow_.clear();
	sasa_blocks_shadow_.resize( pose_size, utility::vector0<uint8_t>(max_rotamer_atoms_) );
	atom_sasa_score_shadow_.clear();
	atom_sasa_score_shadow_.resize( pose_size, utility::vector0<float>(max_rotamer_atoms_) );
	atom_sasa_score_fast_shadow_.clear();
	atom_sasa_score_fast_shadow_.resize( pose_size, &delta_sasa_zeros_.front() );
	atom_sap_shadow_.clear();
	atom_sap_shadow_.resize( pose_size, utility::vector0<float>(max_rotamer_atoms_) );

	current_score_shadow_ = 0;

	shadow_mismatch_.clear();
	shadow_mismatch_.resize( pose_size );
	dirty_sasa_.clear();
	dirty_sasa_.resize( pose_size, utility::vector0<bool>(max_rotamer_atoms_, fast_) );


	lightning_current_res_sap_.clear();
	lightning_current_res_sap_.resize( pose_size );

	lightning_current_res_sap_shadow_.clear();
	lightning_current_res_sap_shadow_.resize( pose_size );

	symm_work_resvect_.clear();
	symm_work_resvect_.resize( pose_size );
	symm_work_resvect_shadow_.clear();
	symm_work_resvect_shadow_.resize( pose_size );

}

void
find_mismatch(
	select::residue_selector::ResidueSubset const & small_sub,
	select::residue_selector::ResidueSubset const & big_sub,
	std::string const & small_name,
	std::string const & big_name
) {
	utility::vector1< Size > small_not_big;
	for ( Size i = 1; i <= small_sub.size(); i++ ) {
		if ( small_sub[i] && ! big_sub[i] ) {
			small_not_big.push_back( i );
		}
	}
	if ( small_not_big.size() > 0 ) {
		TR.Warning << small_name << " selects residues that are not part of " << big_name  << ": " << small_not_big << std::endl;
	}
}



// Resize all arrays to the size of a new pose
void
SapConstraintHelper::resize_arrays(
	core::pose::Pose const & pose
) {

	find_mismatch( score_positions_, sap_calculate_positions_, "score_selector", "sap_calculate_selector" );
	find_mismatch( score_positions_, sasa_positions_, "score_selector", "sasa_selector" );
	find_mismatch( sap_calculate_positions_, sasa_positions_, "sap_calculate_selector", "sasa_selector" );

	check_positions_sap_.clear();
	check_positions_sap_.resize( pose.size() );
	check_positions_block_.clear();
	check_positions_block_.resize( pose.size() );


	delta_sasa_scratch_.clear();
	delta_sasa_scratch_.resize( lightning_ ? 1 : max_rotamer_atoms_);
	delta_sasa_zeros_.clear();
	delta_sasa_zeros_.resize( lightning_ ? 1 : max_rotamer_atoms_);


}

// Look through the SapDatabase for each of the residue types in the current rotamer_sets. The SapDatabase
//  actually contains information on each atom type. So what we do here is we construct each residues type
//  from each atom type in the SapDatabase
//  The block_params then contain enough information to convert from blocks into a sasa_score
void
SapConstraintHelper::fill_block_params( core::pose::Pose const & pose, pack::rotamer_set::RotamerSets const & rotamer_sets ) {
	max_rotamer_atoms_ = 0;
	all_block_params_.clear();
	rotamer_to_block_param_offset_.clear();

	SapDatabase * db = SapDatabase::get_instance();

	// We're going to make the assumtion that if the names are the same the atoms are the same
	std::unordered_map< std::string, Size > res_name_to_offset_;

	for ( Size irot = 1; irot <= rotamer_sets.nrotamers() + pose.size(); irot++ ) {

		conformation::ResidueCOP const & rotamer = irot <= rotamer_sets.nrotamers() ?
			rotamer_sets.rotamer( irot ) :
			pose.residue( irot - rotamer_sets.nrotamers() ).get_self_ptr();

		bool warn = score_positions_[ rotamer->seqpos() ] || sap_calculate_positions_[ rotamer->seqpos() ];
		std::pair< char, std::string > name1_name3 = db->get_name1_name3( *rotamer, warn );
		if ( ! name1_name3.first ) {
			score_positions_[ rotamer->seqpos() ] = false;
			sap_calculate_positions_[ rotamer->seqpos() ] = false;
		}

		// Make sure there are not duplicate rotamers
		debug_assert( rotamer_to_block_param_offset_.count( &*rotamer ) == 0 );

		auto iter = res_name_to_offset_.find( rotamer->name() );

		Size offset = 0;
		if ( iter != res_name_to_offset_.end() ) {
			offset = iter->second;
		} else {

			offset = all_block_params_.size()+1;
			res_name_to_offset_[rotamer->name()] = offset;

			Size first_sidechain = rotamer->first_sidechain_atom();
			Size natoms = my_natoms( rotamer );

			max_rotamer_atoms_ = std::max<Size>( natoms, max_rotamer_atoms_ );

			float total_sasa = 0;
			utility::vector1< BlockParam > block_params;

			// Identify all the atom types and accumulate the max sasa
			for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
				Size iatom = iat + first_sidechain;
				std::string atom_type = name1_name3.second + "_" + utility::trim( rotamer->atom_name(iatom) );

				auto iter2 = db->atomtype_to_block_param()->find( atom_type );
				if ( iter2 == db->atomtype_to_block_param()->end() ) {
					if ( name1_name3.first && warn ) TR.Warning << "Atom type not found: " << atom_type << std::endl;
					block_params.emplace_back( 0, 0, 0 );
				} else {
					block_params.push_back( iter2->second );
					if ( ! rotamer->atom_is_backbone( iatom ) ) {
						total_sasa += iter2->second.max_sasa_score;
					}
				}
			}


			Real hydrophobic_weight = name1_name3.first? my_hydrophobic_weight( db, name1_name3.first ) : 0;

			for ( BlockParam param : block_params ) {
				all_block_params_.emplace_back( param.max_sasa_score / total_sasa * hydrophobic_weight,
					param.no_block, param.full_block );
			}
		}
		rotamer_to_block_param_offset_[ &*rotamer ] = offset;
	}
}



// Calculate the "block" parameter. This is basically just a sigmoid using 3 lines.
//  So a max value, a min value, and a linear interpolation
//  This value is sort of like a pair-wise decomposable sasa
Real
SapConstraintHelper::calculate_block( Real distance, Real radius_us, Real radius_them ) {
	Real full_block = radius_us + radius_them;
	Real full_unblock = full_block + 3;

	if ( distance > full_unblock ) return 0;

	Real block_magnitude = 5 * radius_us*radius_us;

	if ( distance < full_block ) return block_magnitude;

	Real block = (full_unblock - distance) / ( full_unblock - full_block ) * block_magnitude;

	return block;

}

// This function is a doozy, this prepares all of the tables that keep track of nearby atoms
// Most notably, it needs to take care of two tables:
//  1. atom_within_5_ -- this is just a map of which atoms on which rotamers are within 5A of other atoms
//  2. interacting_block_ -- this is a table that gives the rotamer-rotamer pre-calculated blocks for sasa
//
// This function is actually the slow step of the "fast" protocol
void
SapConstraintHelper::fill_atom_neighbor_stuff(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets
) {
	atom_within_5_map_.clear();
	atom_within_5_.clear();
	rotamer_to_sasa_data_.clear();
	interacting_block_offset_.clear();
	interacting_block_.clear();

	const Real R2 = 5*5;
	const Real max_block2 = MAX_BLOCK_INTERACTION*MAX_BLOCK_INTERACTION;



	// If we're in fast, we need to zero out the per-rotamer sasa tables because they aren't standard containers
	// In "slow" we can't pre-calculate this
	if ( fast_ ) {
		for ( Size irot = 1; irot <= rotamer_sets.nrotamers() + pose.size(); irot++ ) {

			conformation::ResidueCOP const & rotamer = irot <= rotamer_sets.nrotamers() ?
				rotamer_sets.rotamer( irot ) :
				pose.residue( irot - rotamer_sets.nrotamers() ).get_self_ptr();

			float * pose_blocks = rotamer_to_sasa_data_[ &*rotamer ];
			std::fill( pose_blocks, pose_blocks + ATOM_SASA_SCORE_ELEMS, 0 );
		}
	}



	pack::task::PackerTaskOP task = rotamer_sets.task()->clone();
	utility::vector1<bool> designing_residues = task->designing_residues();


	// Make a bunch of vectors because this seemed to make things faster when profiling
	utility::vector1<float> fake_blocks( max_rotamer_atoms_ + 10 );

	utility::vector1<Real> iblocks( max_rotamer_atoms_ + 10 );
	utility::vector1<Real> jblocks( max_rotamer_atoms_ + 10 );

	utility::vector1<utility::vector1<uint8_t>> isap( max_rotamer_atoms_ + 10 );
	utility::vector1<utility::vector1<uint8_t>> jsap( max_rotamer_atoms_ + 10 );

	utility::vector1<numeric::xyzVector<Real>> iatoms( max_rotamer_atoms_ + 10 );
	utility::vector1<numeric::xyzVector<Real>> jatoms( max_rotamer_atoms_ + 10 );

	// figure out the max rotamers at a given position so we can size the next two vectors
	Size max_rotamers = 0;
	for ( Size a = 1; a <= pose.size(); a++ ) {
		if ( rotamer_sets.has_rotamer_set_for_residue( a ) ) {
			max_rotamers = std::max<Size>( rotamer_sets.rotamer_set_for_residue( a )->num_rotamers(), max_rotamers );
		}
	}
	utility::vector1<conformation::Residue const *> irotamers( max_rotamers );
	utility::vector1<conformation::Residue const *> jrotamers( max_rotamers );

	// This is the one place we have to desymmetrize the pose because otherwise the packer_graph doesn't work right
	core::pose::Pose desymm_pose;
	if ( symm_info_ ) {
		desymm_pose = pose;
		core::pose::symmetry::make_asymmetric_pose( desymm_pose );
	}

	// This enables design on residue 1. We need to do this so that we can trick create_packer_graph into rebuilding the graph
	task->clean_residue_task( pose.residue(1), 1, pose );

	// This is a convenient way to get all pairs of interacting positions at 2 + 2 + 5 A (the max block interaction distance)
	utility::graph::GraphOP graph = pack::create_packer_graph( symm_info_ ? desymm_pose : pose, *fake_lr_scorefxn_, task );

	// Use this for positions that don't have rotsets
	pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

	// This is the loop where we loop over all positions
	for ( Size i = 1; i <= graph->num_nodes(); i++ ) {
		if ( ! sasa_positions_[i] ) continue;

		// Figure out all positions j that this position i will interact with
		utility::vector1<Size> inspect_positions;
		inspect_positions.push_back( i );
		for ( auto
				iter = graph->get_node( i )->const_upper_edge_list_begin(),
				iter_end = graph->get_node( i )->const_upper_edge_list_end();
				iter != iter_end; ++iter ) {
			Size potential_j = (*iter)->get_other_ind(i);
			if ( ! sasa_positions_[ potential_j ] ) continue;
			inspect_positions.push_back( potential_j );
		}

		// Get the rotamers for this position
		pack::rotamer_set::RotamerSetCOP const & irotset = rotamer_sets.has_rotamer_set_for_residue( i ) ?
			rotamer_sets.rotamer_set_for_residue( i )
			: fake_rotset;

		Size inum_rotamers = irotset->num_rotamers();
		for ( Size a = 1; a <= inum_rotamers; a++ ) {
			irotamers[a] = &*irotset->rotamer( a );
		}


		for ( Size j : inspect_positions ) {

			runtime_assert( i <= j );

			bool positions_had_sap = false;
			bool positions_had_block = false;

			// Get the rotamers for this position
			pack::rotamer_set::RotamerSetCOP const & jrotset = rotamer_sets.has_rotamer_set_for_residue( j ) ?
				rotamer_sets.rotamer_set_for_residue( j )
				: fake_rotset;
			Size jnum_rotamers = jrotset->num_rotamers();
			for ( Size a = 1; a <= jnum_rotamers; a++ ) {
				jrotamers[a] = &*jrotset->rotamer( a );
			}


			// Here we begin the all-rotamer by all-rotamer for loop. We need to calculate which atoms are within 5A of each
			//  other as well as calculate the blocks that each rotamer applies to the other
			for ( Size irot = 0; irot <= inum_rotamers; irot++ ) {

				// Pull out all the information for the i rotamer
				conformation::Residue const * irotamer = irot > 0 ? irotamers[irot] : &pose.residue( i );
				Size ifirst_sidechain = irotamer->first_sidechain_atom();
				Size inatoms = my_natoms( irotamer );
				for ( Size iatom = 1; iatom <= inatoms; iatom++ ) {
					iatoms[iatom] = irotamer->xyz( iatom );
				}
				Distance irad = irotamer->nbr_radius();

				// This is where we store the sasa/block information during "fast" calculations. For "slow" we throw this away
				float * irot_pose_blocks = !fast_ ? &fake_blocks.front() : rotamer_to_sasa_data_.at( irotamer );


				for ( Size jrot = 0; jrot <= jnum_rotamers; jrot++ ) {
					if ( i == j && irot != jrot ) continue; // if it's the same seqpos it has to be the same rotamer

					// Pull out the information for the j rotamer
					conformation::Residue const * jrotamer = jrot > 0 ? jrotamers[jrot] : &pose.residue( j );
					Size jfirst_sidechain = jrotamer->first_sidechain_atom();
					Distance jrad = jrotamer->nbr_radius();

					// Do a quick check to see if we can skip this pair. The block distance is always greater than 5A
					Distance interaction = irad + jrad + MAX_BLOCK_INTERACTION;
					Distance interaction2 = interaction*interaction;
					Distance cb_dist2 = irotamer->nbr_atom_xyz().distance_squared( jrotamer->nbr_atom_xyz() );
					if ( interaction2 < cb_dist2 ) {
						continue;   // so far away that not even the blocks will reach
					}
					// Since in "fast" we don't calculate blocks for rotamers other than 0, we can skip this calculation
					//  if it's clear that no two atoms are within 5A
					if ( fast_ ) {
						if ( irot != 0 && jrot != 0 ) {
							Distance interaction_sap = irad + jrad + 5;
							Distance interaction_sap2 = interaction_sap*interaction_sap;
							if ( interaction_sap2 < cb_dist2 ) {
								continue;
							}
						}
					}

					// This is where we store the sasa/block information during "fast" calculations. For "slow" we throw this away
					float * jrot_pose_blocks = !fast_ ? &fake_blocks.front() : rotamer_to_sasa_data_.at( jrotamer );

					Size jnatoms = my_natoms( jrotamer );
					for ( Size jatom = 1; jatom <= jnatoms; jatom++ ) {
						jatoms[jatom] = jrotamer->xyz( jatom );
					}

					// Clear out our temporary per-atom blocks
					std::fill( iblocks.begin(), iblocks.begin()+inatoms+1, 0 );
					std::fill( jblocks.begin(), jblocks.begin()+jnatoms+1, 0 );

					// Clear out sap atom buffers
					for ( Size a = 1; a <= inatoms; a++ ) {
						isap[a].clear();
					}
					for ( Size a = 1; a <= jnatoms; a++ ) {
						jsap[a].clear();
					}

					// Use these booleans so that we know if we can discard the entire calculation at the end
					bool any_blocks = false;
					bool any_sap = false;

					// This is the true atom-atom loop where we calculate blocks and find pairs within 5A
					// We keep track of backbone atoms because backbone atoms are not part of sap calculations per the paper
					for ( Size iatom = 1; iatom <= inatoms; iatom++ ) {
						bool ibb = irotamer->atom_is_backbone( iatom );

						for ( Size jatom = 1; jatom <= jnatoms; jatom++ ) {
							bool jbb = jrotamer->atom_is_backbone( jatom );

							Real dist2 = iatoms[iatom].distance_squared(jatoms[jatom]);

							// The block radius is larger than 5A so skip early if we can
							if ( dist2 > max_block2 ) continue;

							// "Fast" only does block calculations based off of rotamer 0
							if ( !fast_ || irot == 0 || jrot == 0 ) {

								Real dist = std::sqrt( dist2 );

								Real iblock = calculate_block( dist, irotamer->atom_type(iatom).lj_radius(),
									jrotamer->atom_type(jatom).lj_radius());
								Real jblock = calculate_block( dist, jrotamer->atom_type(jatom).lj_radius(),
									irotamer->atom_type(iatom).lj_radius());

								if ( !fast_ ) {
									// We're storing the entire block of one rotamer to another on a per atom basis
									// This is where we are accumulating it atom by atom
									iblocks[iatom] += iblock;
									jblocks[jatom] += jblock;

									// Backbone atoms don't have sasa so we don't need to remember if they got blocked
									any_blocks |= (iblock > 0 && !ibb) || (jblock > 0 && !jbb);
								} else {
									// For fast, we want to pretend that the original pose is poly CYS if it's designable

									if ( ! designing_residues[i] || (int)iatom - (int)ifirst_sidechain <= 0 ) {
										jblocks[jatom] += jblock;
										any_blocks |= jblock > 0 && !jbb;
									}

									if ( ! designing_residues[j] || (int)jatom - (int)jfirst_sidechain <= 0 ) {
										iblocks[iatom] += iblock;
										any_blocks |= iblock > 0 && !ibb;
									}

								}
							}

							// This is where we store the atom-atom 5A pairs for later
							// This effectively keeps backbone atoms from calculating sap because they won't be added here
							if ( !ibb && !jbb && dist2 < R2 ) {

								isap[iatom].push_back( jatom - jfirst_sidechain );
								jsap[jatom].push_back( iatom - ifirst_sidechain );

								any_sap = true;
							}
						}
					}

					if ( any_blocks && i != j ) {

						if ( ! fast_ ) {

							// Here we fill up interacting_block_ with all sidechain atoms of i and then all sidechain atoms of j
							// This way, we have pre-calculated the complete blockage of each atom by the other residue
							Size offset = interacting_block_.size()+1;
							std::pair< conformation::Residue const *, conformation::Residue const * > key( &*irotamer, &*jrotamer );
							debug_assert( interacting_block_offset_.count( key ) == 0 );
							interacting_block_offset_[ key ] = offset;


							for ( Size iatom = ifirst_sidechain; iatom <= inatoms; iatom++ ) {
								interacting_block_.push_back( (uint8_t)std::min<Size>( 255, lround( iblocks[iatom] / SAP_BLOCK_STORE_SCALE ) ) );
							}
							for ( Size jatom = jfirst_sidechain; jatom <= jnatoms; jatom++ ) {
								interacting_block_.push_back( (uint8_t)std::min<Size>( 255, lround( jblocks[jatom] / SAP_BLOCK_STORE_SCALE ) ) );
							}

						} else {

							// For fast, if this rotamer was calculated against rotamer 0 at the other position, we can store the blocks
							// outright into the sasa table because we only base it off rotamer 0. Later we convert these into sasa scores

							// If these UBs are not on the last atom, then it's a bad rotamer anyways because all canonicals fit
							if ( jrot == 0 ) {
								Size ub = std::min<Size>( inatoms, ifirst_sidechain + ATOM_SASA_SCORE_ELEMS - 1 );
								for ( Size iatom = ifirst_sidechain; iatom <= ub; iatom++ ) {
									irot_pose_blocks[ iatom - ifirst_sidechain ] += iblocks[iatom];
								}
							}
							if ( irot == 0 ) {
								Size ub = std::min<Size>( jnatoms, jfirst_sidechain + ATOM_SASA_SCORE_ELEMS - 1 );
								for ( Size jatom = jfirst_sidechain; jatom <= ub; jatom++ ) {
									jrot_pose_blocks[ jatom - jfirst_sidechain ] += jblocks[jatom];

								}
							}
						}

						positions_had_block |= true;
					}

					// If there was any sap, we store the 5A atom pairs into atom_within_5
					if ( any_sap && sap_calculate_positions_[ i ] && sap_calculate_positions_[ j ] ) {

						add_sap_data( &*irotamer, &*jrotamer, isap, inatoms, ifirst_sidechain );
						add_sap_data( &*jrotamer, &*irotamer, jsap, jnatoms, jfirst_sidechain );

						positions_had_sap = true;
					}


				}
			}

			// These two arrays allow us to quickly know which other seqpos we need to check when looking for block
			//  partners or sap partners
			if ( positions_had_sap ) {
				check_positions_sap_[i].push_back(j);
				if ( i != j ) {
					check_positions_sap_[j].push_back(i);
				}
			}

			if ( positions_had_block ) {
				check_positions_block_[i].push_back(j);
				check_positions_block_[j].push_back(i);
			}
		}
	}

}

// This is the function where we actually store that atom_within_5 data
// The hardest part about this function is the crazy format of atom_within_5
void
SapConstraintHelper::add_sap_data(
	conformation::Residue const * key1,
	conformation::Residue const * key2,
	utility::vector1<utility::vector1<uint8_t>> const & sap_positions,
	Size natoms,
	Size first_sidechain
) {
	// atom_within_5 tells you which atoms on key.second are within 5Å of each atom on key.first
	// Remember: This is the structure of atom_within_5_value
	//           std::pair< Size, uint8_t[ATOM_WITHIN_5_ELEMS]>
	// The uint8 has two sections:
	//   1. A list of length ( natoms - first_sidechain + 1 ) that contains the number of key.second atoms for each atom
	//   2. A bunch of lists of key.second atoms (one list for each atom on key.first but all concated together)
	//
	// There are only 40 bytes to store 2. however. Once it overflows, the Size refers to the offset into
	//   atom_within_5_ where the list continues

	std::pair< conformation::Residue const *, conformation::Residue const * > key( key1, key2 );

	// Just use the dictionary to make it because we can't copy the array
	atom_within_5_value & value = atom_within_5_map_[ key ];

	// Data is the data we're going to store in 2. above. Make it all as 1 big array and then put it into place
	std::vector<uint8_t> data;

	// Length of 1. above
	Size length_entries = natoms - first_sidechain + 1;

	// Fill in 2. above keeping track of how many elements we're adding and storing that in 1.
	for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
		Size iatom = iat + first_sidechain;

		utility::vector1<uint8_t> const & this_data = sap_positions[iatom];
		value.second[iat] = (uint8_t)this_data.size();

		data.insert( data.end(), this_data.begin(), this_data.end() );
	}

	// Actually put the 2. data into place
	if ( data.size() + length_entries < ATOM_WITHIN_5_ELEMS ) {
		// Here the data just fits into atom_within_5_value so we store it there
		value.first = 0;
		for ( Size i = 0; i < data.size(); i++ ) value.second[length_entries+i] = data[i];
	} else {
		// Here the data does not fit and so we first fill in atom_within_5_value as much as we can and then overflow to atom_within_5_
		value.first = atom_within_5_.size()+1;
		for ( Size i = 0; i < ATOM_WITHIN_5_ELEMS-length_entries; i++ ) value.second[length_entries+i] = data[i];
		for ( Size i = ATOM_WITHIN_5_ELEMS-length_entries; i < data.size(); i++ ) atom_within_5_.push_back( data[i] );
	}


}

// Get number of atoms whether we are in heavy_ mode or not
Size
SapConstraintHelper::my_natoms(
	core::conformation::ResidueCOP const & res
) const {
	if ( heavy_ ) {
		return res->nheavyatoms();
	} else {
		return res->natoms();
	}
}

// Get number of atoms whether we are in heavy_ mode or not
Size
SapConstraintHelper::my_natoms(
	core::conformation::Residue const * res
) const {
	if ( heavy_ ) {
		return res->nheavyatoms();
	} else {
		return res->natoms();
	}
}

Real
SapConstraintHelper::my_hydrophobic_weight(
	SapDatabase * db,
	char aa
) const {
	Real weight = db->hydrophobic_weight( aa );
	if ( ! utility::isnan( options_->sap_parameter_options().hydrop_lys_arg_setting ) ) {
		if ( aa == 'R' || aa == 'K' ) {
			weight = options_->sap_parameter_options().hydrop_lys_arg_setting;
		}
	}

	weight += options_->sap_parameter_options().hydrop_adder;

	return weight;
}


// This is a helper function for "fast" since we are filling in the sasa's ahead of time
// All we need to do is step through each atom on each rotamer and use the BlockParams to convert
//  the blocks into a sasa score
void
SapConstraintHelper::convert_block_sasa_to_sasa_score(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets
) {

	for ( Size irot = 1; irot <= rotamer_sets.nrotamers() + pose.size(); irot++ ) {

		conformation::ResidueCOP const & rotamer = irot <= rotamer_sets.nrotamers() ?
			rotamer_sets.rotamer( irot ) :
			pose.residue( irot - rotamer_sets.nrotamers() ).get_self_ptr();


		Size first_sidechain = rotamer->first_sidechain_atom();
		Size natoms = my_natoms( rotamer );

		Size block_param_offset = rotamer_to_block_param_offset_.at( &*rotamer );

		float * sasa_ptr = rotamer_to_sasa_data_[ &*rotamer ];

		int ub = std::min<int>( (int)natoms - (int)first_sidechain, ATOM_SASA_SCORE_ELEMS - 1 );
		for ( Size iat = 0; (int)iat <= ub ; iat++ ) {

			BlockParam const & param = all_block_params_[block_param_offset+iat];

			float sasa;
			uint16_t blocks = sasa_ptr[iat] / SAP_BLOCK_STORE_SCALE;
			float frac_sasa = 0;
			if ( blocks > param.full_block ) {
				sasa = 0;
			} else if ( blocks < param.no_block ) {
				sasa = param.max_sasa_score;
				//frac_sasa = 1;
			} else {
				frac_sasa = float( param.full_block - blocks ) / float( param.full_block - param.no_block );
				sasa = frac_sasa * param.max_sasa_score;
			}

			sasa_ptr[iat] = sasa;
		}

	}

}

// This is a bit of extra work that "fast" needs to do instead of "slow"
// We mark all atoms as dirty permanently so that they always get their saps recalculated
// We also convert the precomputed blocks into precomputed sasa_scores
void
SapConstraintHelper::setup_fast(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets
) {

	for ( Size i = 1; i <= dirty_sasa_.size(); i++ ) {
		for ( Size j = 0; j < max_rotamer_atoms_; j++ ) {
			dirty_sasa_[i][j] = true;
		}
	}

	convert_block_sasa_to_sasa_score( pose, rotamer_sets );
}


// This function performs the sasa_score accumulation from from_rot onto to_rot for lightning
// Since we sum the entire sidechain for sap and since the sasa's are precomputed, we can converge
//  one rotamer's entire interaction to another into one number
Real
SapConstraintHelper::find_lightning_2b(
	conformation::ResidueCOP const & from_rot,
	conformation::ResidueCOP const & to_rot
) {
	Real twobody = 0;

	std::pair< conformation::Residue const *, conformation::Residue const * > key( &*from_rot, &*to_rot );

	Size from_first_sidechain = from_rot->first_sidechain_atom();
	Size from_natoms = my_natoms( from_rot );

	auto iter = atom_within_5_map_.find( key );

	if ( iter == atom_within_5_map_.end() ) return twobody;

	// Get the sasa_scores of the from_rot
	float * from_sasa_ptr = rotamer_to_sasa_data_[ &*from_rot ];

	// atom_within_5 tells you which atoms on key.second are within 5Å of each atom on key.first
	// Remember: This is the structure of atom_within_5_value
	//           std::pair< Size, uint8_t[ATOM_WITHIN_5_ELEMS]>
	// The uint8 has two sections:
	//   1. A list of length ( natoms - first_sidechain + 1 ) that contains the number of key.second atoms for each atom
	//   2. A bunch of lists of key.second atoms (one list for each atom on key.first but all concated together)
	//
	// There are only 40 bytes to store 2. however. Once it overflows, the Size refers to the offset into
	//   atom_within_5_ where the list continues

	// Get atom_within_5_value
	atom_within_5_value const & value = iter->second;

	// offset is offset into atom_within_5_
	// local offset is offset into value.second
	Size offset = value.first;
	Size local_offset = from_natoms - from_first_sidechain + 1;

	Size iat = 0;   // which iat on key.first
	Size elems = 0; // how many atoms on key.second for this iat
	Size ielem = 0; // which atom in elems are we on
	float delta_amount = 0;

	// This would only trigger if a ligand got through, but ligands get flagged as bad and never calculate sap
	debug_assert( (int)from_natoms - (int)from_first_sidechain < (int)ATOM_SASA_SCORE_ELEMS - 1 );

	// Here we iterate inside value.second;
	for ( ; (int)iat <= (int)from_natoms - (int)from_first_sidechain && local_offset < ATOM_WITHIN_5_ELEMS; iat++ ) {

		elems = value.second[iat];
		if ( elems > 0 ) {

			delta_amount = from_sasa_ptr[iat];

			for ( ielem = 0; ielem < elems && local_offset < ATOM_WITHIN_5_ELEMS; ielem++, local_offset++ ) {

				// We give every atom within 5A the sasa score from this atom
				twobody += delta_amount;
			}
		}
	}

	// This is a no-op but we're leaving it for symmetry between sections like this
	offset += local_offset - ATOM_WITHIN_5_ELEMS;

	// clean up the atom we were interrupted on
	for ( ; ielem < elems; ielem++, offset++ ) {

		// We give every atom within 5A the sasa score from this atom
		twobody += delta_amount;
	}

	// Here we iterate inside atom_within_5_
	for ( ; (int)iat <= (int)from_natoms - (int)from_first_sidechain; iat++ ) {

		elems = value.second[iat];
		if ( elems > 0 ) {

			delta_amount = from_sasa_ptr[iat];

			for ( ielem = 0; ielem < elems; ielem++, offset++ ) {

				// We give every atom within 5A the sasa score from this atom
				twobody += delta_amount;
			}
		}
	}
	return twobody;
}

// Lookup into lightning_rotamer_1b_sap_ for a given aa type and seqpos
float &
SapConstraintHelper::lightning_1b_lookup( chemical::AA aa, Size seqpos ) {
	Size aa_idx = lightning_aa_2_index_[ aa ];
	Size seqpos_idx = lightning_seqpos_2_index_[ seqpos ];

	debug_assert( aa_idx != BAD_AA_INDEX );
	debug_assert( seqpos_idx != BAD_SEQPOS_INDEX );

	return lightning_rotamer_1b_sap_[ seqpos_idx * lightning_num_aas_ + aa_idx ];
}



// you'll probably have to draw a picture to understand this one
Size
upper_triangle_offset( Size y, Size x, Size size ) {
	debug_assert( y < x );
	return size * y - ( y * y + y ) / 2 + x - y - 1;
}

// Lookup into lightning_rotamer_2b_sap_ for two given aa types and seqposs
std::pair< float, float > &
SapConstraintHelper::lightning_2b_lookup( chemical::AA aa1, Size seqpos1, chemical::AA aa2, Size seqpos2, bool create /*= false*/ ) {
	debug_assert( seqpos1 < seqpos2 );


	Size aa_idx1 = lightning_aa_2_index_[ (Size)aa1 ];
	Size aa_idx2 = lightning_aa_2_index_[ (Size)aa2 ];


	if ( use_2b_map_ ) {

		basic::datacache::ResRotPair rrp( seqpos1, aa_idx1, seqpos2, aa_idx2 );

		if ( create ) {
			return lightning_rotamer_2b_sap_map_[ rrp ];
		} else {
			auto iter = lightning_rotamer_2b_sap_map_.find( rrp );
			if ( iter == lightning_rotamer_2b_sap_map_.end() ) {
				return default_2b_;
			} else {
				return iter->second;
			}
		}
		runtime_assert( false );
		return default_2b_;

	} else {


		// we only store the upper triangle for seqpos
		// first figure out the offset to the seqpos upper triangle
		// then go inside and figure out the index to the aa

		Size seqpos_idx1 = lightning_seqpos_2_index_[ seqpos1 ];
		Size seqpos_idx2 = lightning_seqpos_2_index_[ seqpos2 ];

		debug_assert( aa_idx1 != BAD_AA_INDEX );
		debug_assert( aa_idx2 != BAD_AA_INDEX );
		debug_assert( seqpos_idx1 != BAD_SEQPOS_INDEX );
		debug_assert( seqpos_idx2 != BAD_SEQPOS_INDEX );

		Size seqpos_part = upper_triangle_offset( seqpos_idx1, seqpos_idx2, lightning_num_seqpos_ );
		Size aa_part = aa_idx1 * lightning_num_aas_ + aa_idx2;

		// std::cout << aa_idx1 << " " << aa_idx2 << " " << seqpos_idx1 << " " << seqpos_idx2 << " " << seqpos_part << " " << aa_part << " " << lightning_rotamer_2b_sap_.size() << " " << lightning_num_seqpos_ << " " << lightning_num_aas_ << std::endl;
		// std::cout << seqpos_part << " " <<  lightning_num_aas_sq_ << " " << aa_part << " " << seqpos_part * lightning_num_aas_sq_ + aa_part << std::endl;

		return lightning_rotamer_2b_sap_[ seqpos_part * lightning_num_aas_sq_ + aa_part ];
	}
}

// Setting up lightning is actually rather similar to setting up the other modes with the exception being that
//  we only use the first rotamer of each aa type at each position
core::pack::rotamer_set::RotamerSets
SapConstraintHelper::setup_lightning(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets
) {
	lightning_rotamer_1b_sap_.clear();
	lightning_rotamer_2b_sap_.clear();
	use_2b_map_ = false;
	lightning_rotamer_2b_sap_map_.clear();
	lightning_aa_2_index_.clear();
	lightning_aa_2_index_.resize( chemical::num_aa_types, BAD_AA_INDEX );
	lightning_num_aas_ = 0;
	lightning_seqpos_2_index_.clear();
	lightning_seqpos_2_index_.resize( pose.size(), BAD_SEQPOS_INDEX );
	lightning_num_seqpos_ = 0;

	// Make rotamer sets containing only 1st rotamer of each aa type at each position

	// Use this for positions that don't have a rotset
	pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

	// lean_rotamer_sets is like rotamer_sets but only contains one rotamer of each aa type at each position
	pack::rotamer_set::RotamerSets lean_rotamer_sets;
	lean_rotamer_sets.set_task( rotamer_sets.task() );

	// Use this set to make sure we only have one of each aa at each position
	std::set< std::pair< Size, Size > > seen_aa_seqpos;

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		pack::rotamer_set::RotamerSetOP new_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();
		new_rotset->set_resid( seqpos );

		// Keep track of which seqpos we are using so we can make the smallest tables possible later
		if ( sap_calculate_positions_[seqpos] ) lightning_seqpos_2_index_[ seqpos ] = lightning_num_seqpos_++;

		// Get the original rotamer set or a blank one
		pack::rotamer_set::RotamerSetCOP old_rotset = rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
			rotamer_sets.rotamer_set_for_residue( seqpos ) : fake_rotset;

		// This for-loop is either a loop over the rotamers or a single loop with just the input pose's residue
		//  It's necessary to use the input pose's residue because if a position isn't packing, it won't have any rotamers
		Size start_rot = 0; //rotamer_sets.has_rotamer_set_for_residue( seqpos ) ? 1 : 0;
		for ( Size irot = start_rot; irot <= old_rotset->num_rotamers(); irot++ ) {
			conformation::ResidueCOP const & rotamer = irot > 0 ? old_rotset->rotamer( irot ) : pose.residue( seqpos ).get_self_ptr();

			// Ensure this is the first rotamer at this position
			std::pair< Size, Size > lookup( (Size)rotamer->aa(), seqpos );
			if ( seen_aa_seqpos.count( lookup ) ) continue;
			seen_aa_seqpos.insert( lookup );

			// Add it to the aa_2_index_ if it hasn't been added yet
			if ( lightning_aa_2_index_.at( (Size)rotamer->aa() ) == BAD_AA_INDEX ) {
				lightning_aa_2_index_[ (Size)rotamer->aa() ] = lightning_num_aas_++;
			}

			new_rotset->add_rotamer( *rotamer->clone() );
		}

		// Add our new rotset to the lean_rotamer_sets
		if ( rotamer_sets.has_rotamer_set_for_residue( seqpos ) ) {
			lean_rotamer_sets.set_explicit_rotamers( lean_rotamer_sets.resid_2_moltenres( seqpos ), new_rotset );
		}
	}

	lean_rotamer_sets.update_offset_data();

	// Now we know lightning_num_aas_ and lightning_num_seqpos_ so we can size these
	lightning_rotamer_1b_sap_.resize( lightning_num_aas_ * lightning_num_seqpos_ );

	lightning_num_aas_sq_ = lightning_num_aas_ * lightning_num_aas_;

	// We'll make 20 GB the flip point
	// 20e9 = num_aa^2 * ( x^2 - x )^2 / 2
	// x ~= 100 when 20 aa
	Size lim = Size(  450 / std::sqrt( lightning_num_aas_ + 1 )  );
	use_2b_map_ = ( lightning_num_seqpos_ > lim ) || symm_debug_force_map_;
	if ( use_2b_map_ ) {
		default_2b_.first = 0;
		default_2b_.second = 0;
		TR << "Using lightning map. Select fewer calculate positions for faster calulations." << std::endl;
	} else {
		lightning_rotamer_2b_sap_.resize( lightning_num_aas_sq_ *
			( lightning_num_seqpos_ * lightning_num_seqpos_ - lightning_num_seqpos_ ) / 2 );
	}


	// This exactly mirrors how we set up "fast", except now we're only doing it with one rotamer of each aa type at each position
	fill_block_params( pose, lean_rotamer_sets );
	fill_atom_neighbor_stuff( pose, lean_rotamer_sets );
	convert_block_sasa_to_sasa_score( pose, lean_rotamer_sets );

	// Calculate pairwise and onebody sap scores

	// Here we fill in the 1b and 2b energies by looping through all pairs of rotamers
	// This is where the sap interactions between each residue are calculated
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sap_calculate_positions_[seqpos] ) continue;

		pack::rotamer_set::RotamerSetOP seq_rotamer_set = lean_rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
			lean_rotamer_sets.rotamer_set_for_residue( seqpos ) : fake_rotset;


		// Loop through all positions that are known to have sap interactions at this position
		utility::vector1< Size > check_positions = check_positions_sap_[ seqpos ];
		for ( Size position : check_positions_sap_[seqpos] ) {
			if ( ! sap_calculate_positions_[position] ) continue;
			if ( position < seqpos ) continue; // don't score backwards so we dont double count

			pack::rotamer_set::RotamerSetOP pos_rotamer_set = lean_rotamer_sets.has_rotamer_set_for_residue( position ) ?
				lean_rotamer_sets.rotamer_set_for_residue( position ) : fake_rotset;

			// Use the same trick from above where we either loop over the rotset or the input pose's residue
			Size start_seq_irot = lean_rotamer_sets.has_rotamer_set_for_residue( seqpos ) ? 1 : 0;
			for ( Size seq_irot = start_seq_irot; seq_irot <= seq_rotamer_set->num_rotamers(); seq_irot++ ) {
				conformation::ResidueCOP const & seq_rot = seq_irot > 0 ? seq_rotamer_set->rotamer( seq_irot ) : pose.residue( seqpos ).get_self_ptr();

				// Same trick
				Size start_pos_irot = lean_rotamer_sets.has_rotamer_set_for_residue( position ) ? 1 : 0;
				for ( Size pos_irot = start_pos_irot; pos_irot <= pos_rotamer_set->num_rotamers(); pos_irot++ ) {
					conformation::ResidueCOP const & pos_rot = pos_irot > 0 ? pos_rotamer_set->rotamer( pos_irot ) : pose.residue( position ).get_self_ptr();

					if ( seqpos == position && pos_rot->aa() != seq_rot->aa() ) continue; // only do the diagonal for self

					// These are the accumulated sasa_scores for each rotamer onto the other
					Real pos_onto_seq = find_lightning_2b( pos_rot, seq_rot );
					Real seq_onto_pos = find_lightning_2b( seq_rot, pos_rot );

					// Store the accumulated saps either into the 1b table or the 2b table
					if ( seqpos == position && pos_rot->aa() == seq_rot->aa() ) {
						lightning_1b_lookup( seq_rot->aa(), seqpos ) = pos_onto_seq;

					} else {
						debug_assert( lightning_2b_lookup( seq_rot->aa(), seqpos, pos_rot->aa(), position ).first == 0 );
						debug_assert( lightning_2b_lookup( seq_rot->aa(), seqpos, pos_rot->aa(), position ).second == 0 );
						// std::cout << "P " << seq_rot->name1() << " " << seqpos  << " " << pos_rot->name1() << " " << position << " -- "  << pos_onto_seq << " " << seq_onto_pos << std::endl;

						if ( use_2b_map_ && pos_onto_seq == 0 && seq_onto_pos == 0 ) continue;

						lightning_2b_lookup( seq_rot->aa(), seqpos, pos_rot->aa(), position, true ) =
							std::pair< float, float>( pos_onto_seq,seq_onto_pos );

					}
				}
			}
		}
	}
	return lean_rotamer_sets;

}



// Called immediately after construction
// This just sets up things that don't change based on input pose
void
SapConstraintHelper::init() {
	// These numbers should all be divisible by 64
	runtime_assert( sizeof(type_atom_within_5_map::value_type) == 64 );
	runtime_assert( sizeof(atom_sasa_score_map::value_type) == 192 );

	fast_ = options_->fast();
	lightning_ = options_->lightning();
	heavy_ = false;
	symm_debug_ = SapDatabase::get_instance()->symm_debug();
	symm_debug_force_map_ = SapDatabase::get_instance()->symm_debug_force_map();

	if ( lightning_ ) runtime_assert( fast_ );

	// We use a fake score function to calculate our interaction graph
	//  we never actually invoke this, we just need it to pass to a function
	fake_lr_scorefxn_ = utility::pointer::make_shared<scoring::ScoreFunction>();
	fake_lr_scorefxn_->set_weight(scoring::fa_rep, 1);
	scoring::methods::EnergyMethodOptions opts = fake_lr_scorefxn_->energy_method_options();
	opts.etable_options().max_dis = MAX_BLOCK_INTERACTION;

	// We initialize all of our tables with this so that we know they aren't initialized
	fake_rotamer_ = core::conformation::get_residue_from_name1('A');

	recalc_positions_scratch_.clear();
	recalc_positions_scratch_.push_back(0);
}

void
SapConstraintHelper::apply_residue_selectors( core::pose::Pose const & pose ) {

	score_positions_ = options_->score_selector()->apply( pose );
	sap_calculate_positions_ = options_->sap_calculate_selector()->apply( pose );
	sasa_positions_ = options_->sasa_selector()->apply( pose );

	// We need to clear out VIRT residues
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( pose.residue( seqpos ).is_virtual_residue() ) {
			score_positions_[ seqpos ] = false;
			sap_calculate_positions_[ seqpos ] = false;
			sasa_positions_[ seqpos ] = false;
		}
	}
}

core::pack::rotamer_set::RotamerSets
SapConstraintHelper::setup_for_symmetry(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotsets
) {
	core::pack::rotamer_set::RotamerSets sym_rotsets;
	symm_info_ = nullptr;
	symm_rotamer_to_other_rotamers_.clear();

	if ( ! core::pose::symmetry::is_symmetric( pose ) ) return sym_rotsets;

	conformation::symmetry::SymmetricConformation const & symm_conf =
		dynamic_cast< conformation::symmetry::SymmetricConformation const & >( pose.conformation() );

	symm_info_ = symm_conf.Symmetry_Info();

	// These subsets only point to the asu and we need to fix that
	utility::vector1< bool > repacking_residues = rotsets.task()->repacking_residues();
	utility::vector1< bool > designing_residues = rotsets.task()->designing_residues();

	// follows_me is a nice shortcut from asu->else
	std::unordered_map< Size, utility::vector1<Size> > follows_me;
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		Size i_follow = symm_info_->bb_follows( seqpos );
		if ( i_follow == 0 ) continue;
		// runtime_assert( i_follow <= seqpos );
		follows_me[ i_follow ].push_back( seqpos );

		// We need to mirror over the packer task stuff too
		repacking_residues[ seqpos ] = repacking_residues[ i_follow ];
		designing_residues[ seqpos ] = designing_residues[ i_follow ];
	}

	// Making packer tasks is terribly verbose
	core::pack::task::operation::OperateOnResidueSubsetOP prevent_task = utility::pointer::make_shared<core::pack::task::operation::OperateOnResidueSubset>(
		utility::pointer::make_shared<core::pack::task::operation::PreventRepackingRLT>(),
		core::select::get_residue_selector_from_subset( repacking_residues ), true );
	core::pack::task::operation::OperateOnResidueSubsetOP restrict_task = utility::pointer::make_shared<core::pack::task::operation::OperateOnResidueSubset>(
		utility::pointer::make_shared<core::pack::task::operation::RestrictToRepackingRLT>(),
		core::select::get_residue_selector_from_subset( designing_residues ), true );

	pack::task::TaskFactoryOP tf = utility::pointer::make_shared< pack::task::TaskFactory >();
	tf->push_back( prevent_task );
	tf->push_back( restrict_task );
	pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );

	sym_rotsets.set_task( task );

	repacking_residues = task->repacking_residues();
	designing_residues = task->designing_residues();

	// We need to call a "static" method on this guy
	core::pack::rotamer_set::symmetry::SymmetricRotamerSet_ static_symm_rotset;

	// If the pose didn't go through the SymmetricRotamerSet_ init routine, it's Tsymm_ is invalid
	//  This is a bad method to fix it.
	core::pose::Pose pose_w_recalc = pose;
	static_cast< core::conformation::symmetry::SymmetricConformation & >( pose_w_recalc.conformation() ).recalculate_transforms();



	bool any_rotsets = false;
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

		Size asu_pos = symm_info_->bb_follows( seqpos );
		if ( asu_pos != 0 ) continue;  // We want to copy from asu to else

		// First, make sure we know how to deal with the input pose itself
		conformation::Residue const & pose_rotamer = pose.residue( seqpos );
		symm_rotamer_to_other_rotamers_[ &pose_rotamer ];   // Trigger creation of the key in case there are no followers
		for ( Size sympos : follows_me[ seqpos ] ) {
			symm_rotamer_to_other_rotamers_[ &pose_rotamer ].push_back( pose.residue( sympos ).get_self_ptr() );
		}


		if ( rotsets.has_rotamer_set_for_residue( seqpos ) ) {
			core::pack::rotamer_set::RotamerSetCOP asu_rotset = rotsets.rotamer_set_for_residue( seqpos );

			// Look, const cast is bad, but I don't really see a way around this.
			//  We need to make a RotamerSets that has the exact same Residues as the input RotamerSets and
			//  if we create a new RotamerSet, it clones all the rotamers on input.
			core::pack::rotamer_set::RotamerSetOP non_const_asu_rotset = std::const_pointer_cast< core::pack::rotamer_set::RotamerSet >( asu_rotset );
			sym_rotsets.set_explicit_rotamers( sym_rotsets.resid_2_moltenres( seqpos ), non_const_asu_rotset );

			// Trigger creation of key in case there are no followers
			for ( Size irot = 1; irot <= asu_rotset->num_rotamers(); irot++ ) symm_rotamer_to_other_rotamers_[ &* asu_rotset->rotamer( irot ) ];

			for ( Size sympos : follows_me[ seqpos ] ) {
				core::pack::rotamer_set::RotamerSetOP sym_rotset = utility::pointer::make_shared<core::pack::rotamer_set::RotamerSet_>();
				sym_rotset->set_resid( sympos );

				for ( Size irot = 1; irot <= asu_rotset->num_rotamers(); irot++ ) {
					conformation::ResidueCOP const & rotamer = asu_rotset->rotamer( irot );
					conformation::ResidueOP new_rotamer = static_symm_rotset.orient_rotamer_to_symmetric_partner( pose_w_recalc, *rotamer, sympos, true );
					new_rotamer->seqpos( sympos );
					sym_rotset->add_rotamer( *new_rotamer );
					conformation::ResidueCOP new_rotamer_cloned = sym_rotset->rotamer( irot ); // It got cloned so we just pull back the cloned version
					symm_rotamer_to_other_rotamers_[ &*rotamer ].push_back( new_rotamer_cloned );
				}
				sym_rotsets.set_explicit_rotamers( sym_rotsets.resid_2_moltenres( sympos ), sym_rotset );

				any_rotsets = true;
			}
		} else {
			// std::cout << "No rotamer set " << seqpos << std::endl;
		}
	}

	if ( any_rotsets ) {
		sym_rotsets.update_offset_data();
	}

	// now we fix the selections

	// first clear out the symmetric copies
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( symm_info_->bb_follows( seqpos ) == 0 ) continue;
		sasa_positions_[ seqpos ] = false;
		sap_calculate_positions_[ seqpos ] = false;
		score_positions_[ seqpos ] = false;
	}

	// Now mirror the asu only for calculate
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( symm_info_->bb_follows( seqpos ) != 0 ) continue;
		for ( Size sympos : follows_me.at( seqpos ) ) {
			sasa_positions_[ sympos ] = sasa_positions_[ seqpos ];
			sap_calculate_positions_[ sympos ] = sap_calculate_positions_[ seqpos ];
			if ( symm_debug_ ) {
				score_positions_[ sympos ] = score_positions_[ seqpos ];
			}
		}
	}

	if ( symm_debug_ ) TR << "symm_debug_ activated" << std::endl;

	return sym_rotsets;
}


// Set up all the stuff we need to perform the calculations from a fresh pose
core::pack::rotamer_set::RotamerSets
SapConstraintHelper::init_with_pose(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets_in
) {
	apply_residue_selectors( pose );
	core::pack::rotamer_set::RotamerSets sym_rotsets = setup_for_symmetry( pose, rotamer_sets_in );
	core::pack::rotamer_set::RotamerSets const & rotamer_sets = symm_info_ ? sym_rotsets : rotamer_sets_in;

	if ( ! lightning_ ) {
		fill_block_params( pose, rotamer_sets );
	}
	resize_arrays( pose );

	if ( ! lightning_ ) {
		fill_atom_neighbor_stuff( pose, rotamer_sets );

		if ( fast_ ) {
			setup_fast( pose, rotamer_sets );
		}
	} else {
		setup_lightning( pose, rotamer_sets );
	}

	return rotamer_sets;
}



SapConstraintOptionsCOP
SapConstraintHelper::options() const {
	return options_;
}


void
SapConstraintHelper::report_final_score( Real actual_sap ) const {

	TR << boost::str(boost::format("Internal sap: %.1f Actual sap: %.1f Suggested packing_correction: %.1f")%current_score_shadow_
		%actual_sap%(actual_sap-current_score_shadow_)) << std::endl;;
}


Real
SapConstraintHelper::current_score() const {
	return current_score_shadow_;
}


//////////////////////////////////////////////////////////////////////////////
/////////////////// HACKY FUNCTIONS /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// This is used when we want to calculate sap correctly. All of the functions in this class
//  work like the real sap calculation except for the sasa part. Here we use an actual sasa
//  calculation to set the sasas and then calculate the real sap score
//  The helper is in a totally invalid state after this call
Real
SapConstraintHelper::set_accurate_sasa_and_recalc( pose::Pose const & pose ) {

	// You have to call calculate_energy() first and this ensures it was called
	runtime_assert( pose.size() == internal_resvect_.size() );
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sap_calculate_positions_[seqpos] ) continue;
		// not a perfect test but will work 90% of the time
		// can't check pointers because the resvect is cloned in rotamer sets
		runtime_assert( pose.residue(seqpos).name3() == internal_resvect_[seqpos]->name3() );
	}


	core::id::AtomID_Map<Real> atom_sasa = sap_atom_sasa( pose, sasa_positions_ );

	SapDatabase * db = SapDatabase::get_instance();

	utility::vector1< Size > recalc_positions;
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sap_calculate_positions_[seqpos] ) continue;

		conformation::Residue const & res = pose.residue(seqpos);

		std::pair< char, std::string > name1_name3 = db->get_name1_name3( res, false );

		Size first_sidechain = res.first_sidechain_atom();
		Size natoms = my_natoms( res.get_self_ptr() );

		Real max_sasa = name1_name3.first ? db->max_sasa( name1_name3.first ) : 0;
		Real hydro = name1_name3.first ? my_hydrophobic_weight( db, name1_name3.first ) : 0;

		for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
			Size iatom = iat + first_sidechain;

			Real score = atom_sasa( seqpos, iatom ) / max_sasa * hydro;
			atom_sasa_score_[seqpos][iat] = score;
			atom_sap_[seqpos][iat] = 0;
		}
		recalc_positions.push_back( seqpos );
	}
	current_score_ = 0;

	// Since this function already invalidates the helper, I don't feel too bad doing this
	// In order to allow all speeds to calculate sap correctly, we have to disable fast_
	// for this function call. fast and lightning both have the required pre-calculations to
	// calculate sap correctly if you tell them to (and set the sasa scores)
	bool old_fast = fast_;
	fast_ = false;
	recalculate_saps( recalc_positions );
	fast_ = old_fast;
	save_to_shadow();

	return current_score_;
}


core::id::AtomID_Map<Real>
SapConstraintHelper::get_per_atom_sap( pose::Pose const & pose ) const {

	// You have to call calculate_energy() first and this ensures it was called
	runtime_assert( pose.size() == internal_resvect_.size() );
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sap_calculate_positions_[seqpos] ) continue;
		// not a perfect test but will work 90% of the time
		// can't check pointers because the resvect is cloned in rotamer sets
		runtime_assert( pose.residue(seqpos).name3() == internal_resvect_[seqpos]->name3() );
	}

	core::id::AtomID_Map<Real> saps;
	core::pose::initialize_atomid_map( saps, pose, Real(0) );

	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sap_calculate_positions_[seqpos] ) continue;
		if ( ! score_positions_[seqpos] ) continue;

		core::conformation::Residue const & res = pose.residue(seqpos);
		Size first_sidechain = res.first_sidechain_atom();
		Size natoms = my_natoms( &res );
		for ( int iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
			Size iatom = iat + first_sidechain;
			saps( seqpos, iatom ) = atom_sap_[seqpos][iat];
		}
	}
	return saps;
}

void
SapConstraintHelper::report() {
	restore_from_shadow();
	for ( Size i = 1; i <= internal_resvect_.size(); i++ ) {
		if ( ! sap_calculate_positions_[i] ) continue;
		Real sum = 0;
		core::conformation::ResidueCOP const & rotamer = internal_resvect_[i];
		Size first_sc = rotamer->first_sidechain_atom();
		for ( Size j = first_sc; j <= my_natoms(rotamer); j++ ) {
			Real value = atom_sap_[i][j-first_sc];
			if ( value > 0 ) sum += value;
		}
		std::cout << i << internal_resvect_[i]->name3() << " " << sum << " -";
		for ( Size j = first_sc; j <= my_natoms(rotamer); j++ ) {
			Real value = atom_sap_[i][j-first_sc];
			std::cout << " " << value << std::endl;
		}
		std::cout << std::endl;
	}
}




} //sap
} //guidance_scoreterms
} //pack
} //core
