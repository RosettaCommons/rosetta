// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

	virtual
	void build_rotamers(
		pose::Pose const & the_pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	);


	/// @brief  Build rotamers that depend on positions of rotamers built in a previous pass
	virtual
	void build_dependent_rotamers(
		RotamerSets const & rotamer_sets,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph
	);


	virtual
	void
	add_rotamer(
		conformation::Residue const & rotamer
	);

	virtual
	Size
	get_n_residue_types() const;

	virtual
	Size
	get_n_residue_groups() const;

	virtual
	Size
	get_residue_type_begin( Size which_restype ) const;

	virtual
	Size
	get_residue_group_begin( Size which_resgroup ) const;

	virtual
	Size
	get_n_rotamers_for_residue_type( Size which_resgroup ) const;

	virtual
	Size
	get_n_rotamers_for_residue_group( Size which_resgroup ) const;

	virtual
	Size
	get_residue_type_index_for_rotamer( Size which_rotamer ) const ;

	virtual
	Size
	get_residue_group_index_for_rotamer( Size which_rotamer ) const ;


	/// @brief Computes the packers "one body energies" for the set of rotamers.
	virtual
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & energies ) const;


	/// for OPTE
	virtual
	void
	compute_one_body_energy_maps(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		utility::vector1< scoring::EnergyMap > & energies ) const;

	/*
	virtual
	void
	compute_two_body_energies(
		RotamerSet const & other,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		ObjexxFCL::FArray2< Energy > & pair_energy_table ) const;
	*/

	virtual
	Size
	num_rotamers() const;

	virtual
	Size
	id_for_current_rotamer() const;

	virtual
	conformation::ResidueCOP
	rotamer( Size rot_id ) const;

	virtual Rotamers::const_iterator begin() const { return rotamers_.begin(); }
	virtual Rotamers::const_iterator end() const { return rotamers_.end(); }

	virtual
	conformation::ResidueOP
	nonconst_rotamer( Size rot_id );

	virtual
	void
	store_trie( Size method_enum_id, conformation::AbstractRotamerTrieOP trie );

	virtual
	conformation::AbstractRotamerTrieCOP
	get_trie( Size method_enum_id ) const;

	/// @brief removes a single rotamer and causes a rotamer index update
	virtual
	void
	drop_rotamer( Size rot_id );

	/// @brief rotamers_to_delete must be of size nrotmaers -- each position
	/// in the array that's "true" is removed from the set of rotamers
	virtual
	void
	drop_rotamers( utility::vector1< bool > const & rotamers_to_delete );

	/// @brief deletes the rotamers in the list with the given indices.
	/// The indices of these rotamers is presumed to be those before any delete operation.
	/// e.g. if there are four rotamers, and rotamer_indices_to_delete includes 1 & 3,
	/// then the rotamers that will remain are the rotamers originally indexed as 2 and 4,
	/// even though their new indices will be 1 & 2.
	virtual
	void
	drop_rotamers_by_index( utility::vector1< Size > const & rotamer_indices_to_delete );

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
		graph::GraphCOP packer_neighbor_graph,
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
		graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
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
		graph::GraphCOP packer_neighbor_graph
	);


	/// @brief Creates a sets of rotamers for an "optimize H" repacking
	void
	build_optimize_H_rotamers(
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		graph::GraphCOP packer_neighbor_graph
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
		graph::GraphCOP packer_neighbor_graph
	) const;

private:

	/// @brief logic for building TP3 water rotamers
	void build_tp3_water_rotamers(
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		graph::GraphCOP packer_neighbor_graph
	);

	/// @brief declare that a new block of residue types has begun, and that new residues
	/// are about to be pushed back.
	void
	prepare_for_new_residue_type( core::chemical::ResidueType const & restype );

	/// @brief should two residue types be considered the same residue type?
	bool
	different_restype( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const;

	/// @brief should two residue types be considered to belong to the same residue-type group?
	bool
	different_resgroup( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const;

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

	void
	update_rotamer_offsets() const;

// DATA
private:

	BumpSelector bump_selector_;
	Rotamers rotamers_;

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
};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

