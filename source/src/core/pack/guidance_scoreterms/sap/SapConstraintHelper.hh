// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh
/// @brief  Score term that applies the RotamerPSSMConstraint as an energy
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintHelper_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintHelper_hh

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapDatabase.hh>

// Package headers
#include <core/id/AtomID_Map.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <basic/datacache/CacheableResRotPairFloatMap.hh>

// Utility headers
#include <utility/vector1.hh>

#include <unordered_map>

// forward declaration for testing
class SapConstraintEnergyTests;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

template<typename Tval>
struct MyTemplatePointerHash1 {
	size_t operator()(const Tval* val) const {
		static const size_t shift = (size_t)log2(1 + sizeof(Tval));
		return (size_t)(val) >> shift;
	}
};


struct RotamerRotamerHasher {
	// std::hash<conformation::Residue const *> hasher;
	MyTemplatePointerHash1<conformation::Residue > hasher;

	std::size_t operator()( std::pair< conformation::Residue const *, conformation::Residue const * > const & k ) const {

		return (hasher( k.first ) << 1) ^ hasher( k.second );
	}
};


// struct StringSizeHasher {
//     // std::hash<conformation::Residue const *> hasher;
//     std::hash<std::string> string_hasher;

//     std::size_t operator()( std::pair< std::string, Size > const & k ) const {

//         return string_hasher( k.first ) ^ k.second;
//     }
// };

// // This is symmetric with respect to first and second because I know I'll never have them switched
// struct StringSizeStringSizeHasher {
//     // std::hash<conformation::Residue const *> hasher;
//     StringSizeHasher string_size_hasher;

//     std::size_t operator()( std::pair< std::pair< std::string, Size >, std::pair< std::string, Size > > const & k ) const {

//         return string_size_hasher( k.first ) ^ string_size_hasher( k.second );
//     }
// };

#define ATOM_WITHIN_5_ELEMS 40
#define ATOM_SASA_SCORE_ELEMS 46

#define BAD_AA_INDEX 65535
#define BAD_SEQPOS_INDEX 4294967295

class SapConstraintHelper {
	friend class ::SapConstraintEnergyTests;

public:

	SapConstraintHelper( SapConstraintOptionsCOP const & options );

	core::Real
	calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const &resvect,
		core::Size const substitution_position
	);

	void
	commit_considered_substitution();


	SapConstraintOptionsCOP
	options() const;

	core::pack::rotamer_set::RotamerSets
	init_with_pose(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	Real current_score() const;

	Real
	set_accurate_sasa_and_recalc( pose::Pose const & pose );

	void report();

	void report_final_score( Real actual_sap ) const;

	core::id::AtomID_Map<Real>
	get_per_atom_sap( pose::Pose const & pose ) const;

protected:

	typedef std::pair< Size, uint8_t[ATOM_WITHIN_5_ELEMS]> atom_within_5_value;
	typedef std::unordered_map< std::pair< conformation::Residue const *, conformation::Residue const * >,
		atom_within_5_value, RotamerRotamerHasher > type_atom_within_5_map;

	typedef float atom_sasa_score_value[ATOM_SASA_SCORE_ELEMS] ;
	typedef std::unordered_map< conformation::Residue const*, atom_sasa_score_value > atom_sasa_score_map;


	void
	restore_from_shadow();

	void
	save_to_shadow();

	core::Real
	symm_calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const &resvect,
		core::Size const substitution_position
	);


	void
	recalculate_saps( utility::vector1<Size> const & positions_to_update );

	void
	recalculate_sasa();

	void
	update_neighbors_sap(
		float const * delta_sasa_score,
		Size seqpos,
		bool invert,
		bool skip_self
	);

	void
	add_to_sap( float amount, Size seqpos, Size iat );

	void
	add_remove_rotamer(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		Size const substitution_position,
		bool add
	);

	void
	add_remove_rotamer_fast(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		Size const substitution_position,
		bool add
	);

	void
	add_remove_rotamer_lightning(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		Size const substitution_position,
		bool add
	);

	void
	lightning_update_neighbors_sap(
		Size seqpos,
		bool invert,
		bool dont_update_self
	);

	void
	lightning_add_to_sap(
		float amount,
		Size seqpos
	);

	void store_sasa_blocks( Size & offset, Size seqpos, bool add, bool mark_dirty );

	void
	reinit_with_resvect(
		utility::vector1< core::conformation::ResidueCOP > const &resvect,
		bool skip_reset = false
	);

	void
	reset_calculation();


	void
	apply_residue_selectors( core::pose::Pose const & pose );

	void
	resize_arrays(
		core::pose::Pose const & pose
	);

	void
	fill_block_params(
		core::pose::Pose const & pose,
		pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	Real
	calculate_block( Real dist2, Real radius_us, Real radius_them );

	void
	fill_atom_neighbor_stuff(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	void
	convert_block_sasa_to_sasa_score(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	float &
	lightning_1b_lookup( chemical::AA aa, Size seqpos );

	std::pair< float, float > &
	lightning_2b_lookup( chemical::AA aa1, Size seqpos1, chemical::AA aa2, Size seqpos2, bool create = false );

	Real
	find_lightning_2b(
		conformation::ResidueCOP const & from_rot,
		conformation::ResidueCOP const & to_rot
	);

	core::pack::rotamer_set::RotamerSets
	setup_lightning(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	void
	setup_fast(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	);

	void
	add_sap_data(
		conformation::Residue const * key1,
		conformation::Residue const * key2,
		utility::vector1<utility::vector1<uint8_t>> const & sap_positions,
		Size natoms,
		Size first_sidechain
	);

	core::pack::rotamer_set::RotamerSets
	setup_for_symmetry(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotsets
	);

	Size
	my_natoms(
		core::conformation::ResidueCOP const & res
	) const;

	Size
	my_natoms(
		core::conformation::Residue const * res
	) const;


	Real
	my_hydrophobic_weight(
		SapDatabase * db,
		char aa
	) const;


	void
	init();


	// void
	// assert_saps_total_score( std::string const & tag);


private:


	//// Description of variables

	SapConstraintOptionsCOP options_;

	// Set the mode of the whole calculation. See top of SapConstraintHelper.cc for explanation
	// Heavy isn't used currently
	bool fast_;
	bool lightning_;
	bool heavy_;

	// Residue subsets
	utility::vector1<bool> score_positions_;
	utility::vector1<bool> sap_calculate_positions_;
	utility::vector1<bool> sasa_positions_;


	//// Calculated once

	core::conformation::ResidueCOP fake_rotamer_;

	core::scoring::ScoreFunctionOP fake_lr_scorefxn_;
	core::scoring::ScoreFunctionOP fake_sr_scorefxn_;



	//// Calculated once per pose

	Size max_rotamer_atoms_;

	// For each positions, which other positions are close enough to check
	utility::vector1<utility::vector1<Size>> check_positions_sap_;
	utility::vector1<utility::vector1<Size>> check_positions_block_;

	// List of BlockParam for each rotamer type. Length of each internal vector
	//  is natoms() - first_sidechain_atom() + 1
	// These block params are different from the ones in the database
	// These ones have already been divided by max sasa and multiplied by hydrophobic_weight
	std::unordered_map< conformation::Residue const*, Size > rotamer_to_block_param_offset_;
	utility::vector1<BlockParam> all_block_params_;



	// Used for the "fast" calculation. This is a dictionary of arrays that store the precomputed
	//  sasa scores for each rotamer
	atom_sasa_score_map rotamer_to_sasa_data_;


	// For each rotamer pair, store interacting block scores
	// First a list of blocks on the key.first
	// Then a list of blocks on key.second
	std::unordered_map< std::pair< conformation::Residue const *, conformation::Residue const *>, Size,
		RotamerRotamerHasher > interacting_block_offset_;
	utility::vector1<uint8_t> interacting_block_;


	// Store all atoms within 5 A of each atom for each residue pair
	// Order of residue pair matters here
	// Value.first represents how far into atom_within_5_ to look for the atom information
	// Value.second is the length of the atom information for each atom followed by
	//     the start of the information stored in atom_within_5_
	//
	// We already pay the cache-line penalty for the dict lookup so we store as much information there as possible
	type_atom_within_5_map atom_within_5_map_;
	utility::vector1<uint8_t> atom_within_5_;



	//// Things that have shadows

	// Copy of the current residues that we have
	utility::vector1< core::conformation::ResidueCOP > internal_resvect_;

	// Offset to all_block_params_ for the current rotamer at each position
	utility::vector1<Size> block_param_offset_;

	// The next few data structures are always valid except while being updated. In that case
	// sasa_blocks_ are updated followed by atom_sasa_score_ at the same positions
	//  Then for each of the atom_sasa_score_ that changes, all atom_sap_ that depend
	//  on that position are updated too (which also updates current_score_)
	// For the most part += and -= are always used

	// Lowest level of data storage. Updated whenever rotamers are added and removed.
	// These are accumulated into atom_sasa_score_ during recalculate_sasa()
	//  dirty_sasa_ exists to indicate that sasa_blocks_ and atom_sasa_score_ are out of synce
	utility::vector1<utility::vector0<uint16_t>> sasa_blocks_;

	// Second lowest level of storage. These are calculated directly from sasa_blocks_
	// These are calculated in recalculate_sasa() and the delta value is used to update atom_sap_
	utility::vector1<utility::vector0<float>> atom_sasa_score_;

	// Third lowest level of storage. These are almost always updated during recalculate_sasa()
	//  with the exception being when a rotamer is first added, they must be fully recalculated.
	utility::vector1<utility::vector0<float>> atom_sap_;

	// Final level of storage. Updated whenever atom_sap_ changes.
	float current_score_;


	utility::vector1<float *> atom_sasa_score_fast_;

	//// State trackers

	// All the positions that need to be changed during a shadow operation
	utility::vector1<bool> shadow_mismatch_;

	// All positions where sasa_blocks_ has been updated but atom_sasa_score_ hasn't been recalculated.
	utility::vector1<utility::vector0<bool>> dirty_sasa_;


	// The shadows
	utility::vector1<utility::vector0<float>> atom_sasa_score_shadow_;
	utility::vector1<utility::vector0<uint16_t>> sasa_blocks_shadow_;
	utility::vector1<utility::vector0<float>> atom_sap_shadow_;
	float current_score_shadow_;
	utility::vector1<Size> block_param_offset_shadow_;
	utility::vector1< core::conformation::ResidueCOP > internal_resvect_shadow_;

	utility::vector1<float *> atom_sasa_score_fast_shadow_;


	// Stuff for lightning

	// Fast mapping of aa to internal index
	utility::vector1< uint16_t > lightning_aa_2_index_;
	Size lightning_num_aas_;
	Size lightning_num_aas_sq_;

	// Fast mapping of seqpos to internal index
	utility::vector1< uint32_t > lightning_seqpos_2_index_;
	Size lightning_num_seqpos_;

	// The sap score a given aa-seqpos casts upon itself
	utility::vector0< float > lightning_rotamer_1b_sap_;

	// The sap score a given aa-seqpos casts upon another
	// Stored in upper triangle mode
	utility::vector0< std::pair< float, float > > lightning_rotamer_2b_sap_;

	bool use_2b_map_;
	std::unordered_map< basic::datacache::ResRotPair, std::pair< float, float >, basic::datacache::ResRotPairHasher > lightning_rotamer_2b_sap_map_;
	std::pair< float, float > default_2b_;

	// The current sap score of each residue
	utility::vector1< Real > lightning_current_res_sap_;

	utility::vector1< Real > lightning_current_res_sap_shadow_;





	//// Scratch
	utility::vector1<Size> recalc_positions_scratch_;
	utility::vector0<float> delta_sasa_scratch_;

	utility::vector0<float> delta_sasa_zeros_;



	//// Everything that had to be added because of symmetry
	conformation::symmetry::SymmetryInfoCOP symm_info_;
	// Size symm_max_following_;
	std::unordered_map< core::conformation::Residue const *, utility::vector1< core::conformation::ResidueCOP > > symm_rotamer_to_other_rotamers_;
	utility::vector1< core::conformation::ResidueCOP > symm_work_resvect_;
	utility::vector1< core::conformation::ResidueCOP > symm_work_resvect_shadow_;
	bool symm_debug_;
	bool symm_debug_force_map_;

};


} //sap
} //guidance_scoreterms
} //pack
} //core

#endif
