// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/LongRangeEnergyContainer.hh
/// @brief  A container interface for storing and scoring long range energies
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_LREnergyContainer_hh
#define INCLUDED_core_scoring_LREnergyContainer_hh

// Unit headers
#include <core/scoring/LREnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/EnergyMap.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace scoring {

class ResidueNeighborIterator : public utility::pointer::ReferenceCount
{
public:
	~ResidueNeighborIterator() override;

	virtual ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & ) = 0;
	virtual ResidueNeighborIterator const & operator ++ () = 0;
	virtual bool operator == ( ResidueNeighborIterator const & ) const = 0;
	virtual bool operator != ( ResidueNeighborIterator const & ) const = 0;

	virtual Size upper_neighbor_id() const = 0;
	virtual Size lower_neighbor_id() const = 0;

	virtual Size residue_iterated_on() const = 0;
	virtual Size neighbor_id() const = 0;

	virtual void save_energy( EnergyMap const & ) = 0;

	virtual void retrieve_energy( EnergyMap & ) const = 0;
	virtual void accumulate_energy( EnergyMap & ) const= 0;

	virtual void mark_energy_computed() = 0;
	virtual void mark_energy_uncomputed() = 0;

	virtual bool energy_computed() const = 0;
};

class ResidueNeighborConstIterator : public utility::pointer::ReferenceCount
{
public:
	~ResidueNeighborConstIterator() override;

	virtual ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & ) = 0;
	virtual ResidueNeighborConstIterator const & operator ++ () = 0;
	virtual bool operator == ( ResidueNeighborConstIterator const & ) const = 0;
	virtual bool operator != ( ResidueNeighborConstIterator const & ) const = 0;

	virtual Size upper_neighbor_id() const = 0;
	virtual Size lower_neighbor_id() const = 0;

	virtual Size residue_iterated_on() const = 0;
	virtual Size neighbor_id() const = 0;

	virtual void retrieve_energy( EnergyMap & ) const = 0;
	virtual void accumulate_energy( EnergyMap & ) const= 0;

	virtual bool energy_computed() const = 0;

};

class LREnergyContainer : public utility::pointer::ReferenceCount
{
public:

	~LREnergyContainer() override;

	virtual
	LREnergyContainerOP clone() const = 0;

	virtual
	bool empty() const = 0;

	virtual
	void
	set_num_nodes( Size ) { }

	virtual
	bool
	any_neighbors_for_residue( int resid ) const = 0;

	virtual
	bool
	any_upper_neighbors_for_residue( int resid ) const = 0;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const = 0;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const = 0;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const = 0;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const = 0;

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) = 0;

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) = 0;

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid ) = 0;

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) = 0;

};

} // namespace scoring
} // namespace core

#endif
