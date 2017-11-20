// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.hh
/// @brief  rotamer set implementation class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSet__hh
#define INCLUDED_core_pack_rotamer_set_RotamerSet__hh

//Unit headers
#include <core/pack/rotamer_set/RotamerSet_.fwd.hh>

//Package headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/BumpSelector.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

//Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#ifdef WIN32
#include <core/scoring/trie/RotamerTrieBase.hh>
#endif
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh> // wym

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

/// @brief Container for a set of rotamers for use in packing.
/// Rotamers are sorted into groups of the same residue type.
/// Offsets into these rotamer groups are maintained by this class, as is
/// information concerning the "original rotamer" -- the rotamer
/// present on the input pose before packing began.
class RotamerSet_ : public RotamerSet
{
public:
	typedef conformation::ResidueOP ResidueOP;
	typedef conformation::ResidueCOP ResidueCOP;
	typedef pack::dunbrack::ChiSetOP ChiSetOP;
	typedef scoring::trie::RotamerTrieBaseOP RotamerTrieBaseOP;

public:
	RotamerSet_();
	virtual ~RotamerSet_();

	void build_rotamers(
		pose::Pose const & the_pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	) override;


	/// @brief  Build rotamers that depend on positions of rotamers built in a previous pass
	void build_dependent_rotamers(
		RotamerSets const & rotamer_sets,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph
	) override;

	void
	add_rotamer(
		conformation::Residue const & rotamer
	) override;

	void
	add_rotamer_into_existing_group(
		conformation::Residue const & rotamer
	) override;

	Size
	get_n_residue_types() const override;

	Size
	get_n_residue_groups() const override;

	Size
	get_residue_type_begin( Size which_restype ) const override;

	Size
	get_residue_group_begin( Size which_resgroup ) const override;

	Size
	get_n_rotamers_for_residue_type( Size which_resgroup ) const override;

	Size
	get_n_rotamers_for_residue_group( Size which_resgroup ) const override;

	Size
	get_residue_type_index_for_rotamer( Size which_rotamer ) const  override;

	Size
	get_residue_group_index_for_rotamer( Size which_rotamer ) const override;


	/// @brief Computes the packers "one body energies" for the set of rotamers.
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & energies ) const override;

	/// @brief Computes the packers one body energies for the set of rotamers as well
	/// as two body energies for neighboring positions defined as packable by the task.
	void
	compute_one_and_two_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & one_body_energies,
		utility::vector1< utility::vector1< core::PackerEnergy > > & two_body_energies,
		utility::vector1< core::Size > & packable_neighbors ) const override;

	/// for OPTE
	void
	compute_one_body_energy_maps(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< scoring::EnergyMap > & energies ) const override;

	/*
	void
	compute_two_body_energies(
	RotamerSet const & other,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	ObjexxFCL::FArray2< Energy > & pair_energy_table ) const;
	*/

	Size
	num_rotamers() const override;

	Size
	id_for_current_rotamer() const override;

	conformation::ResidueCOP
	rotamer( Size rot_id ) const override;

	basic::datacache::BasicDataCache &
	rotamer_data_cache( Size rot_id ) const override;

	conformation::Residue const &
	rotamer_ref( Size rot_id ) const override;

	Rotamers::const_iterator begin() const override { return rotamers_.begin(); }
	Rotamers::const_iterator end() const override { return rotamers_.end(); }

	conformation::ResidueOP
	nonconst_rotamer( Size rot_id ) override;

	void
	store_trie( Size method_enum_id, conformation::AbstractRotamerTrieOP trie ) override;

	conformation::AbstractRotamerTrieCOP
	get_trie( Size method_enum_id ) const override;

	/// @brief removes a single rotamer and causes a rotamer index update
	void
	drop_rotamer( Size rot_id ) override;

	/// @brief rotamers_to_delete must be of size nrotmaers -- each position
	/// in the array that's "true" is removed from the set of rotamers
	void
	drop_rotamers( utility::vector1< bool > const & rotamers_to_delete ) override;

	/// @brief deletes the rotamers in the list with the given indices.
	/// The indices of these rotamers is presumed to be those before any delete operation.
	/// e.g. if there are four rotamers, and rotamer_indices_to_delete includes 1 & 3,
	/// then the rotamers that will remain are the rotamers originally indexed as 2 and 4,
	/// even though their new indices will be 1 & 2.
	void
	drop_rotamers_by_index( utility::vector1< Size > const & rotamer_indices_to_delete ) override;

	/// @brief Give the pose a chance to stash any data needed by the _rotset_
	///        need nonconst access to pose
	void
	initialize_pose_for_rotset_creation(
		pose::Pose & /*pose*/
	) const override {}

private:
	RotamerSet_( RotamerSet_ const & );


protected:
	/// @brief Creates a set of rotamers for a particular residue type
	/// (the concrete residue type) while relying on the rotamer-
	/// building instructions within the PackerTask.
	/// Use the residue in the input pose at position resid_ as the existing residue.
	virtual
	void build_rotamers_for_concrete_virt(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		utility::graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	);

	/// @brief Creates a set of rotamers for a particular residue type
	/// (the concrete residue type) while relying on the rotamer-
	/// building instructions within the PackerTask.
	void build_rotamers_for_concrete(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	);


	/// @brief, Filter water rotamers by eliminating rotamers with overlaps
	/// and uniformly only keeping a maximum number of rotamers
	/// (water_rotamers_cap, def=500)
	void
	filter_water_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< conformation::ResidueOP > const & new_rotamers,
		utility::vector1< conformation::ResidueOP > & filtered_rotamers
	);


	/// @brief  Build rotamers of a specific type that depend on positions of rotamers built in a previous pass
	/// Use an input "existing residue" which may or may not reflect the coordinates of the residue at
	/// pose.residue( resid_ );
	void
	build_dependent_rotamers_for_concrete(
		RotamerSets const & rotamer_sets,
		pose::Pose const & pose,
		scoring::ScoreFunction const &, // scorefxn,
		task::PackerTask const & task,
		conformation::Residue const & existing_residue,
		chemical::ResidueTypeCOP concrete_residue,
		utility::graph::GraphCOP packer_neighbor_graph
	);


	/// @brief Creates a sets of rotamers for an "optimize H" repacking
	void
	build_optimize_H_rotamers(
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::graph::GraphCOP packer_neighbor_graph,
		scoring::ScoreFunction const & scorefxn
	);


public:
	/// @brief Pushes standard-deviation multiples that should be sampled
	/// for this residue -- if this residue has more neighbors within 10A
	/// than the task-specified cutoff for buriedness, then extra rotamer
	/// samples are added to the extra_chi_steps vector, otherwise, the
	/// vector is not modified.
	void set_extra_samples(
		task::PackerTask const & task,
		int num_10A_neighbors,
		int chi,
		chemical::ResidueTypeCOP concrete_residue,
		utility::vector1< Real > & extra_chi_steps
	) const;

public:
	/// @brief Computes the "bump energy" of a rotamer: the bump energy is the
	/// sum of rotamer's interactions with 1) the backbone-and-side chains of
	/// neighboring residues that are held fixed during this repacking optimization
	/// and 2) the backbones of neighboring residues that are changable during this
	/// repacking optimization.
	virtual
	core::PackerEnergy
	bump_check(
		ResidueCOP rotamer,
		scoring::ScoreFunction const & sf,
		pose::Pose const & pose,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph
	) const;

	void
	show( std::ostream & out ) const override;

private:

	/// @brief logic for building TP3 water rotamers
	void build_tp3_water_rotamers(
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::graph::GraphCOP packer_neighbor_graph,
		scoring::ScoreFunction const & scorefxn
	);

	void build_filtered_tp3_water_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::graph::GraphCOP packer_neighbor_graph
	);

	/// @brief declare that a new block of residue types has begun, and that new residues
	/// are about to be pushed back.
	void
	prepare_for_new_residue_type( core::chemical::ResidueType const & restype );

	/// @brief This function should not be called directly -- it ought to be called only from prepare_for_new_residue_type
	void
	new_residue_type();

	/// @brief This function should not be called directly -- it ought to be called only from prepare_for_new_residue_type
	void
	new_residue_group();

	/// @brief appends a rotamer to the list of rotamers, and increments the count
	/// for the number of rotamers for the current value of n_residue_types.
	void
	push_back_rotamer( conformation::ResidueOP );

	/// @brief append a vector of rotamers to the list of rotamers,
	/// and increments the count for the number of rotamers for the current value of n_residue_types.
	/// It is assumed that all of the rotamers in the vector of are the current type and group.
	void
	push_back_rotamers( Rotamers const & );

	/// @brief Lazy update of rotamer indices and offsets and integration of those rotamers
	/// in the rotamers_waiting_for_sort_ list.
	void
	update_rotamer_offsets() const;

	// DATA
private:

	BumpSelector bump_selector_;

	mutable Rotamers rotamers_;
	mutable std::list< ResidueOP > rotamers_waiting_for_sort_;

	mutable Size n_residue_types_;
	mutable Size n_residue_groups_;
	mutable utility::vector1< Size > residue_type_rotamers_begin_;
	mutable utility::vector1< Size > residue_group_rotamers_begin_;
	mutable utility::vector1< Size > n_rotamers_for_restype_;
	mutable utility::vector1< Size > n_rotamers_for_resgroup_;

	mutable utility::vector1< Size > residue_type_for_rotamers_;
	mutable utility::vector1< Size > residue_group_for_rotamers_;

	utility::vector1< conformation::AbstractRotamerTrieOP > cached_tries_;

	mutable Size id_for_current_rotamer_;

	ResidueOP current_rotamer_copy_;
	mutable bool rotamer_offsets_require_update_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief should two residue types be considered the same residue type?
bool
different_restype( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 );

/// @brief should two residue types be considered to belong to the same residue-type group?
bool
different_resgroup( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 );

void
sort_new_rotamers_into_rotset_vector(
	utility::vector1< conformation::ResidueOP > & rotamers,
	std::list< conformation::ResidueOP > & rotamers_waiting_for_sort,
	core::Size & id_for_current_rotamer
);


} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSet_ )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

