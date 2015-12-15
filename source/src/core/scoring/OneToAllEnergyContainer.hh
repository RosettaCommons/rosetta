// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/OneToAllEnergyContainer.hh
/// @brief  A container interface for storing and scoring long range energies
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_OneToAllEnergyContainer_hh
#define INCLUDED_core_scoring_OneToAllEnergyContainer_hh

// Unit headers
#include <core/scoring/OneToAllEnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class OneToAllNeighborIterator : public ResidueNeighborIterator
{
public:
	virtual ~OneToAllNeighborIterator();

	OneToAllNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > * table_in,
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
	Size pos1_;
	Size pos2_;
	bool operating_on_pos1_;
	ScoreType score_type_;
	utility::vector1< Real > * table_;
	utility::vector1< bool > * computed_;
};


///////////////////////////////////////////////////////

class OneToAllNeighborConstIterator : public ResidueNeighborConstIterator
{
public:
	virtual ~OneToAllNeighborConstIterator();

	OneToAllNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > const * table_in,
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
	Size pos1_;
	Size pos2_;
	bool operating_on_pos1_;
	ScoreType score_type_;
	utility::vector1< Real > const * table_;
	utility::vector1< bool > const * computed_;

};

///////////////////////////////////////////////////////////////////////////

class OneToAllEnergyContainer : public LREnergyContainer
{
public:
	virtual
	~OneToAllEnergyContainer();

	virtual
	LREnergyContainerOP clone() const;

	OneToAllEnergyContainer( int const fixed_res_idx, Size const size_in, ScoreType const score_type_in );

	virtual
	bool empty() const;

	virtual
	void
	set_num_nodes( Size size_in );

	virtual
	bool
	any_neighbors_for_residue( int ) const;

	virtual
	bool
	any_upper_neighbors_for_residue( int resid ) const;

	Size
	size() const;

	int
	fixed() const;

	//////////////////// const versions
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
	int fixed_;
	Size size_;
	ScoreType score_type_;

	utility::vector1< Real > table_;
	utility::vector1< bool > computed_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	OneToAllEnergyContainer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_OneToAllEnergyContainer )
#endif // SERIALIZATION


#endif
