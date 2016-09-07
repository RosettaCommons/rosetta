// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DenseEnergyContainer.hh
/// @brief  A container for storing all-against-all energies as the upper triangle of a matrix
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_DenseEnergyContainer_hh
#define INCLUDED_core_scoring_DenseEnergyContainer_hh

// Unit headers
#include <core/scoring/DenseEnergyContainer.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class DenseNeighborIterator : public ResidueNeighborIterator
{
	DenseNeighborIterator & operator = (DenseNeighborIterator const & src );
public:
	~DenseNeighborIterator() override;

	DenseNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > * table_in,
		ObjexxFCL::FArray2D< bool > * computed_in
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
	ScoreType score_type_;
	ObjexxFCL::FArray2D< Real > * table_;
	ObjexxFCL::FArray2D< bool > * computed_;

};


///////////////////////////////////////////////////////

class DenseNeighborConstIterator : public ResidueNeighborConstIterator
{
	DenseNeighborConstIterator & operator = (DenseNeighborConstIterator const & src );
public:
	~DenseNeighborConstIterator() override;

	DenseNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > const * table_in,
		ObjexxFCL::FArray2D< bool > const * computed_in
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
	ScoreType score_type_;
	ObjexxFCL::FArray2D< Real > const * table_;
	ObjexxFCL::FArray2D< bool > const * computed_;

};

///////////////////////////////////////////////////////////////////////////

class DenseEnergyContainer : public LREnergyContainer
{
public:
	~DenseEnergyContainer() override;

	
	LREnergyContainerOP clone() const override;

	DenseEnergyContainer( Size const size_in, ScoreType const score_type_in );

	
	bool empty() const override;

	
	void
	set_num_nodes( Size size_in ) override;

	
	bool
	any_neighbors_for_residue( int ) const override;

	
	bool
	any_upper_neighbors_for_residue( int resid ) const override;

	Size
	size() const;

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
	Size /*const*/ size_;
	ScoreType score_type_;

	ObjexxFCL::FArray2D< Real > table_;
	ObjexxFCL::FArray2D< bool > computed_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DenseEnergyContainer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_DenseEnergyContainer )
#endif // SERIALIZATION


#endif
