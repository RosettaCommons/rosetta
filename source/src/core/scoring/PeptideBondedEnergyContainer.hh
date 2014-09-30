// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PeptideBondedEnergyContainer.hh
/// @brief  A container interface long range energies for n->n+1 interactions only
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_PeptideBondedEnergyContainer_hh
#define INCLUDED_core_scoring_PeptideBondedEnergyContainer_hh

// Unit headers
#include <core/scoring/PeptideBondedEnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class PeptideBondedNeighborIterator : public ResidueNeighborIterator
{
public:
	virtual ~PeptideBondedNeighborIterator();

	PeptideBondedNeighborIterator(
		Size const base_in,
		Size const pos_in,
		utility::vector1< ScoreType > const st,
		utility::vector1< utility::vector1< Real > > * table_in,
		utility::vector1< bool > * computed_in
	);

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & src );

	virtual ResidueNeighborIterator const & operator ++ ();

	virtual bool operator == ( ResidueNeighborIterator const & other ) const;

	virtual bool operator != ( ResidueNeighborIterator const & other ) const;

	virtual Size upper_neighbor_id() const;

	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;

	virtual Size neighbor_id() const;

	virtual void save_energy( EnergyMap const & emap );

	virtual void retrieve_energy( EnergyMap & emap ) const;

	virtual void accumulate_energy( EnergyMap & emap ) const;

	virtual void mark_energy_computed();

	virtual void mark_energy_uncomputed();

	virtual bool energy_computed() const;

private:
	Size base_;
	Size pos_;
	utility::vector1< ScoreType > score_types_;
	utility::vector1< utility::vector1< Real > > * tables_;
	utility::vector1< bool > * computed_;
};


///////////////////////////////////////////////////////

class PeptideBondedNeighborConstIterator : public ResidueNeighborConstIterator
{
public:
	virtual ~PeptideBondedNeighborConstIterator();

	PeptideBondedNeighborConstIterator(
		Size const base_in,
		Size const pos_in,
		utility::vector1< ScoreType > const st,
		utility::vector1< utility::vector1< Real > > const * table_in,
		utility::vector1< bool > const * computed_in
	);

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & src );

	virtual ResidueNeighborConstIterator const & operator ++ ();

	virtual bool operator == ( ResidueNeighborConstIterator const & other ) const;

	virtual bool operator != ( ResidueNeighborConstIterator const & other ) const;

	virtual Size upper_neighbor_id() const;

	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;

	virtual Size neighbor_id() const;

	virtual void retrieve_energy( EnergyMap & emap ) const;

	virtual void accumulate_energy( EnergyMap & emap ) const;

	virtual bool energy_computed() const;

private:
	Size base_;
	Size pos_;
	utility::vector1< ScoreType > score_types_;
	utility::vector1< utility::vector1< Real > > const * tables_;
	utility::vector1< bool > const * computed_;
};

///////////////////////////////////////////////////////////////////////////

class PeptideBondedEnergyContainer : public LREnergyContainer {
public:
	virtual
	~PeptideBondedEnergyContainer();

	virtual
	LREnergyContainerOP clone() const;

	PeptideBondedEnergyContainer( Size const size_in, utility::vector1< ScoreType > const score_type_in, Size offset_in=0 );

	virtual
	bool empty() const;

	using LREnergyContainer::set_num_nodes;

	virtual
	void
	set_num_nodes( Size size_in, Size offset_in=0 );

	virtual
	bool
	any_neighbors_for_residue( int /*resid*/ ) const;

	virtual
	bool
	any_upper_neighbors_for_residue( int /*resid*/ ) const;

	Size
	size() const;

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

	//////////////////// non-const versions
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

private:
	Size size_, offset_;
	utility::vector1< ScoreType > score_types_;
	utility::vector1< utility::vector1< Real > > tables_;
	utility::vector1< bool > computed_;

};

} // namespace scoring
} // namespace core

#endif
