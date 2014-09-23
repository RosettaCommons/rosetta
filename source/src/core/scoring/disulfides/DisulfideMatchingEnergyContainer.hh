// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingEnergyContainer.hh
/// @brief  Disulfide Energy Container class declaration
/// @author rvernon@u.washington.edu
/// @todo This file was taken with minimal modifications from FullatomDisulfideEnergyContainer.hh.
/// It might be a good idea to move some code into a superclass inheritted by both.

#ifndef INCLUDED_core_scoring_disulfides_DisulfideMatchingEnergyContainer_hh
#define INCLUDED_core_scoring_disulfides_DisulfideMatchingEnergyContainer_hh

// Unit headers
#include <core/scoring/disulfides/DisulfideMatchingEnergyContainer.fwd.hh>

#ifdef WIN32
#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#endif

// Package headers
#include <core/scoring/LREnergyContainer.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>

// AUTO-REMOVED #include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

/** @brief An iterator over the disulfide bonds a residue forms
 *
 * When scoring a pose, a long range energy container must be able to iterate
 * over all the residues which interact with a particular residue. For
 * disulfide bonds, this is either zero or one items depending on whether the
 * residue specified forms a disulfide bond or not.
 *
 * @todo Given the proper options, include all residues in the vicinity which
 * could form a bond, not just the best one. Maybe even include non-cysteines!
 */
class DisulfideMatchingNeighborIterator : public ResidueNeighborIterator {
public:

	DisulfideMatchingNeighborIterator(
		DisulfideMatchingEnergyContainer * owner,
		Size focused_node,
		Size disulfide_index
	);

	DisulfideMatchingNeighborIterator(
		DisulfideMatchingEnergyContainer * owner
	);

	virtual ~DisulfideMatchingNeighborIterator();

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & );
	virtual ResidueNeighborIterator const & operator ++ ();
	virtual bool operator == ( ResidueNeighborIterator const & ) const;
	virtual bool operator != ( ResidueNeighborIterator const & ) const;

	virtual Size upper_neighbor_id() const;
	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;
	virtual Size neighbor_id() const;

	virtual void save_energy( EnergyMap const & );
	virtual void retrieve_energy( EnergyMap & ) const;
	virtual void accumulate_energy( EnergyMap & ) const;

	virtual void mark_energy_computed();
	virtual void mark_energy_uncomputed();

	virtual bool energy_computed() const;

private:
	DisulfideMatchingEnergyContainer * owner_;
	Size focused_residue_;
	Size disulfide_index_;
};

/// @brief Just a const version of DisulfideMatchingNeighborIterator
class DisulfideMatchingNeighborConstIterator : public ResidueNeighborConstIterator {
public:
	DisulfideMatchingNeighborConstIterator(
		DisulfideMatchingEnergyContainer const * owner,
		Size focused_node,
		Size disulfide_index
	);

	DisulfideMatchingNeighborConstIterator( DisulfideMatchingEnergyContainer const * owner );

	virtual ~DisulfideMatchingNeighborConstIterator();

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & );
	virtual ResidueNeighborConstIterator const & operator ++ ();
	virtual bool operator == ( ResidueNeighborConstIterator const & ) const;
	virtual bool operator != ( ResidueNeighborConstIterator const & ) const;

	virtual Size upper_neighbor_id() const;
	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;
	virtual Size neighbor_id() const;

	virtual void retrieve_energy( EnergyMap & ) const;
	virtual void accumulate_energy( EnergyMap & ) const;

	virtual bool energy_computed() const;

private:
	DisulfideMatchingEnergyContainer const * owner_;
	Size focused_residue_;
	Size disulfide_index_;

};

/**
 * @brief Storage for Disulfide Energy Terms
 */
class DisulfideMatchingEnergyComponents
{
public:
	DisulfideMatchingEnergyComponents() :
		dslfc_rot_( 0.0 ),
		dslfc_trans_( 0.0 ),
		dslfc_RT_( 0.0 )
	{}

	Energy dslfc_rot() const { return dslfc_rot_;}
	Energy dslfc_trans() const { return dslfc_trans_;}
	Energy dslfc_RT() const { return dslfc_RT_;}

	Energy & dslfc_rot() { return dslfc_rot_;}
	Energy & dslfc_trans() { return dslfc_trans_;}
	Energy & dslfc_RT() { return dslfc_RT_;}

private:
	Energy dslfc_rot_;
	Energy dslfc_trans_;
	Energy dslfc_RT_;
};

class DisulfideMatchingEnergyContainer : public LREnergyContainer {
public:
	static Size const NO_DISULFIDE;

public:
	DisulfideMatchingEnergyContainer();

	DisulfideMatchingEnergyContainer( pose::Pose const & );

	void
	update( pose::Pose const & );

	virtual
	~DisulfideMatchingEnergyContainer();

	virtual
	bool empty() const;

	virtual
	LREnergyContainerOP clone() const;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const;

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid );

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid );

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid );

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid );

	bool
	disulfide_bonded( Size res1id, Size res2id ) const;

	bool residue_forms_disulfide( Size resid ) const;

	// What residue is resid forming a disulfide with?
	Size other_neighbor_id( Size resid ) const;

	DisulfideAtomIndices const &
	disulfide_atom_indices( Size resid ) const;

	DisulfideAtomIndices const &
	other_neighbor_atom_indices( Size resid ) const;

// Read and write access granted to Disulf Iterators. Other classes should not use these methods.
public:

	// Mutators
	void save_energy( Size disulfide_index, EnergyMap const & emap );
	void mark_energy_computed( Size disulfide_index );
	void mark_energy_uncomputed( Size disulfide_index );

	// Accessors
	Size lower_neighbor_id( Size disulfide_index ) const;
	Size upper_neighbor_id( Size disulfide_index ) const;
	Size other_neighbor_id( Size disulfide_index, Size resid ) const;

	void accumulate_energy( Size disulfide_index, EnergyMap & emap ) const;
	void retrieve_energy( Size disulfide_index, EnergyMap & emap ) const;
	bool energy_computed( Size disulfide_index ) const;

private:
	void find_disulfides( pose::Pose const & pose );
	bool disulfides_changed( pose::Pose const & pose );
	Size num_disulfides() const;
private:
	utility::vector1< Size > resid_2_disulfide_index_;
	utility::vector1< chemical::ResidueTypeCOP > disulfide_residue_types_;
	utility::vector1< std::pair< Size, Size > > disulfide_partners_;
	utility::vector1< std::pair< DisulfideAtomIndices, DisulfideAtomIndices > > disulfide_atom_indices_;
	utility::vector1< std::pair< DisulfideMatchingEnergyComponents, bool > > disulfide_info_;
};

}
}
}

#endif
