// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	virtual ~DenseNeighborIterator();

	DenseNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > * table_in,
		ObjexxFCL::FArray2D< bool > * computed_in
	);

	virtual ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & src );

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
	ScoreType score_type_;
	ObjexxFCL::FArray2D< Real > * table_;
	ObjexxFCL::FArray2D< bool > * computed_;

};


///////////////////////////////////////////////////////

class DenseNeighborConstIterator : public ResidueNeighborConstIterator
{
	DenseNeighborConstIterator & operator = (DenseNeighborConstIterator const & src );
public:
	virtual ~DenseNeighborConstIterator();

	DenseNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > const * table_in,
		ObjexxFCL::FArray2D< bool > const * computed_in
	);

	virtual ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & src );

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
	ScoreType score_type_;
	ObjexxFCL::FArray2D< Real > const * table_;
	ObjexxFCL::FArray2D< bool > const * computed_;

};

///////////////////////////////////////////////////////////////////////////

class DenseEnergyContainer : public LREnergyContainer
{
public:
	virtual ~DenseEnergyContainer();

	virtual
	LREnergyContainerOP clone() const;

	DenseEnergyContainer( Size const size_in, ScoreType const score_type_in );

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
