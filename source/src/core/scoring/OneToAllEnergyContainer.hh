// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
	OneToAllNeighborIterator & operator = (OneToAllNeighborIterator const & src );
public:
	~OneToAllNeighborIterator() override;

	OneToAllNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > * table_in,
		utility::vector1< bool > * computed_in
	);

	ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & src ) override;

	ResidueNeighborIterator const & operator ++ () override;

	bool operator == ( ResidueNeighborIterator const & other ) const override;

	bool operator != ( ResidueNeighborIterator const & other ) const override;

	Size upper_neighbor_id() const override;

	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;

	Size neighbor_id() const override;

	void save_energy( EnergyMap const & emap ) override;

	void retrieve_energy( EnergyMap & emap ) const override;

	void accumulate_energy( EnergyMap & emap ) const override;

	void mark_energy_computed() override;

	void mark_energy_uncomputed() override;

	bool energy_computed() const override;

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
	OneToAllNeighborConstIterator & operator = (OneToAllNeighborConstIterator const & src );
public:
	~OneToAllNeighborConstIterator() override;

	OneToAllNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > const * table_in,
		utility::vector1< bool > const * computed_in
	);

	ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & src ) override;

	ResidueNeighborConstIterator const & operator ++ () override;

	bool operator == ( ResidueNeighborConstIterator const & other ) const override;

	bool operator != ( ResidueNeighborConstIterator const & other ) const override;

	Size upper_neighbor_id() const override;

	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;

	Size neighbor_id() const override;

	void retrieve_energy( EnergyMap & emap ) const override;

	void accumulate_energy( EnergyMap & emap ) const override;

	bool energy_computed() const override;

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

	~OneToAllEnergyContainer() override;


	LREnergyContainerOP clone() const override;

	OneToAllEnergyContainer( int const fixed_res_idx, Size const size_in, ScoreType const score_type_in );


	bool empty() const override;


	void
	set_num_nodes( Size size_in ) override;


	bool
	any_neighbors_for_residue( int ) const override;


	bool
	any_upper_neighbors_for_residue( int resid ) const override;

	Size
	size() const;

	int
	fixed() const;

	//////////////////// const versions

	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const override;

	//////////////////// non-const versions

	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) override;


	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) override;


	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid ) override;


	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) override;

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
