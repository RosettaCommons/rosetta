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


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSubset_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSubset_hh

//Unit headers
#include <core/pack/rotamer_set/RotamerSubset.fwd.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>


//Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/graph/Graph.fwd.hh>


// Utility headers
#include <utility/pointer/owning_ptr.hh>

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
	#include <core/pack/rotamer_set/RotamerSet.hh>
	#include <core/graph/Graph.hh>
#endif

namespace core {
namespace pack {
namespace rotamer_set {

/// @brief Container for a subset of rotamers that have been created
/// by another rotamer set.  This subset object copies pointers to
/// the rotamers contained in another set, as opposed to cloning the
/// rotamers.  It's main purpose is to manage the bookkeeping involved
/// in packing with a subset of rotamers (when it might be faster
/// to use a subset and to create an interaction graph specifically
/// for that subset than to simply pass an abreviated list of rotamers
/// to the SimAnnealer with the "rot_to_pack" vector).
class RotamerSubset : public RotamerSet
{
public:
	typedef conformation::ResidueOP ResidueOP;
	typedef conformation::ResidueCOP ResidueCOP;
	typedef scoring::trie::RotamerTrieBaseOP RotamerTrieBaseOP;

public:
	RotamerSubset(
		RotamerSet & rotset,
		utility::vector1< Size > const & rotamer_subset
	);

	virtual ~RotamerSubset();

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
	get_n_rotamers_for_residue_type( Size which_restype ) const;

	virtual
	Size
	get_n_rotamers_for_residue_group( Size which_resgroup ) const;

	/// @brief given a rotamer id, return an int which represents a type for this rotamer.
	virtual
	Size
	get_residue_type_index_for_rotamer( Size which_rotamer ) const ;

	virtual
	Size
	get_residue_group_index_for_rotamer( Size which_rotamer ) const;

	virtual
	Size
	num_rotamers() const;

	virtual
	Size
	id_for_current_rotamer() const;

	virtual
	conformation::ResidueCOP
	rotamer( Size rot_id ) const;

	virtual
	conformation::Residue const &
	rotamer_ref( Size rot_id ) const;

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
	/// @brief (private) No copy-constructor
	RotamerSubset( RotamerSubset const & );

	/// @brief declare that a new block of residue types has begun, and that new residues
	/// are about to be pushed back. NOT IMPLEMENTED.
	//void
	//declare_new_residue_type();

	/// @brief appends a rotamer to the list of rotamers, and increments the count
	/// for the number of rotamers for the current value of n_residue_types. NOT IMPLEMENTED.
	//void
	//push_back_rotamer( conformation::ResidueOP );

	/// @brief Borrow (steal) a rotamer held by another RotamerSet
	/// without cloning that rotamer.  That is, both sets will now point
	/// at the same rotamer object, so if that rotamer changes for one set,
	/// it changes for both.
	void
	steal_rotamer( conformation::ResidueOP rotamer );


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


public: // noop functions:

	virtual
	void build_rotamers(
		pose::Pose const & the_pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	);

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
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & energies ) const;

	virtual
	void
	compute_one_and_two_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,		
		utility::vector1< core::PackerEnergy > & one_body_energies,
		utility::vector1< utility::vector1< core::PackerEnergy > > & two_body_energies,
		utility::vector1< core::Size > & packable_neighbors ) const;
		
	/// for OptE
	virtual
	void
	compute_one_body_energy_maps(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		utility::vector1< scoring::EnergyMap > & energies ) const;


// DATA
private:

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

