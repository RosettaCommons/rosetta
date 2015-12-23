// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/CentroidDisulfideEnergyContainer.hh
/// @brief  Disulfide Energy Container class declaration
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @todo This file was taken with minimal modifications from FullatomDisulfideEnergyContainer.hh.
/// It might be a good idea to move some code into a superclass inheritted by both.

#ifndef INCLUDED_core_scoring_disulfides_CentroidDisulfideEnergyContainer_hh
#define INCLUDED_core_scoring_disulfides_CentroidDisulfideEnergyContainer_hh

// Unit headers
#include <core/scoring/disulfides/CentroidDisulfideEnergyContainer.fwd.hh>

#ifdef WIN32
#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#endif

// Package headers
#include <core/scoring/LREnergyContainer.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>


#include <core/pose/Pose.fwd.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


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
class CentroidDisulfideNeighborIterator : public ResidueNeighborIterator {
public:

	CentroidDisulfideNeighborIterator(
		CentroidDisulfideEnergyContainer * owner,
		Size focused_node,
		Size disulfide_index
	);

	CentroidDisulfideNeighborIterator(
		CentroidDisulfideEnergyContainer * owner
	);

	virtual ~CentroidDisulfideNeighborIterator();

	virtual ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & );
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
	CentroidDisulfideEnergyContainer * owner_;
	Size focused_residue_;
	Size disulfide_index_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CentroidDisulfideNeighborIterator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< CentroidDisulfideNeighborIterator > & construct );
#endif // SERIALIZATION

};

/// @brief Just a const version of CentroidDisulfideNeighborIterator
class CentroidDisulfideNeighborConstIterator : public ResidueNeighborConstIterator {
public:
	CentroidDisulfideNeighborConstIterator(
		CentroidDisulfideEnergyContainer const * owner,
		Size focused_node,
		Size disulfide_index
	);

	CentroidDisulfideNeighborConstIterator( CentroidDisulfideEnergyContainer const * owner );

	virtual ~CentroidDisulfideNeighborConstIterator();

	virtual ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & );
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
	CentroidDisulfideEnergyContainer const * owner_;
	Size focused_residue_;
	Size disulfide_index_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CentroidDisulfideNeighborConstIterator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< CentroidDisulfideNeighborConstIterator > & construct );
#endif // SERIALIZATION

};

/**
* @brief Storage for Disulfide Energy Terms
* @note Although the full atom and centroid terms will not be defined
*       at the same time for a particular disulfide bond, it is convienient
*       to be able to store either in the same object.
*/
class CentroidDisulfideEnergyComponents
{
public:
	CentroidDisulfideEnergyComponents() :
		dslfc_cen_dst_( 0.0 ),
		dslfc_cb_dst_( 0.0 ),
		dslfc_ang_( 0.0 ),
		dslfc_cb_dih_( 0.0 ),
		dslfc_bb_dih_( 0.0 )
	{}

	Energy dslfc_cen_dst() const { return dslfc_cen_dst_;}
	Energy dslfc_cb_dst() const { return dslfc_cb_dst_;}
	Energy dslfc_ang() const { return dslfc_ang_;}
	Energy dslfc_cb_dih() const { return dslfc_cb_dih_;}
	Energy dslfc_bb_dih() const { return dslfc_bb_dih_;}

	Energy & dslfc_cen_dst() { return dslfc_cen_dst_;}
	Energy & dslfc_cb_dst() { return dslfc_cb_dst_;}
	Energy & dslfc_ang() { return dslfc_ang_;}
	Energy & dslfc_cb_dih() { return dslfc_cb_dih_;}
	Energy & dslfc_bb_dih() { return dslfc_bb_dih_;}

private:
	Energy dslfc_cen_dst_;
	Energy dslfc_cb_dst_;
	Energy dslfc_ang_;
	Energy dslfc_cb_dih_;
	Energy dslfc_bb_dih_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class CentroidDisulfideEnergyContainer : public LREnergyContainer {
public:
	static Size const NO_DISULFIDE;

public:
	CentroidDisulfideEnergyContainer();

	CentroidDisulfideEnergyContainer( pose::Pose const & );

	void
	update( pose::Pose const & );

	virtual
	~CentroidDisulfideEnergyContainer();

	virtual
	bool empty() const;

	virtual
	LREnergyContainerOP clone() const;

	virtual
	bool
	any_neighbors_for_residue( int resid ) const;

	virtual
	bool
	any_upper_neighbors_for_residue( int resid ) const;

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
	utility::vector1< std::pair< CentroidDisulfideEnergyComponents, bool > > disulfide_info_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_disulfides_CentroidDisulfideEnergyContainer )
#endif // SERIALIZATION


#endif
