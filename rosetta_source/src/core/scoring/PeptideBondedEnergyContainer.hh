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

namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class PeptideBondedNeighborIterator : public ResidueNeighborIterator
{
public:
	virtual ~PeptideBondedNeighborIterator(){}

	PeptideBondedNeighborIterator(
		Size const pos_in,
		ScoreType const st,
		utility::vector1< Real > * table_in,
		utility::vector1< bool > * computed_in
	):
		pos_( pos_in ),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & src ) {
		assert( dynamic_cast< PeptideBondedNeighborIterator const * >( &src ) );
		PeptideBondedNeighborIterator const & my_src( static_cast< PeptideBondedNeighborIterator const & >( src ) );
		pos_ = my_src.pos_;
		table_ = my_src.table_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborIterator const & operator ++ () {
		++pos_;
		return *this;
	}

	virtual bool operator == ( ResidueNeighborIterator const & other ) const
	{
		return ( residue_iterated_on() == other.residue_iterated_on() );
	}

	virtual bool operator != ( ResidueNeighborIterator const & other ) const
	{
		return !( *this == other );
	}

	virtual Size upper_neighbor_id() const {
		return pos_+1;
	}

	virtual Size lower_neighbor_id() const {
		return pos_;
	}

	virtual Size residue_iterated_on() const {
		return pos_;
	}

	virtual Size neighbor_id() const {
		return pos_+1;
	}

	virtual void save_energy( EnergyMap const & emap ) {
		Real const energy( emap[ score_type_ ] );
		(*table_)[ pos_ ] = energy;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const {
		emap[ score_type_ ] = (*table_)[pos_];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const {
		emap[ score_type_ ] += (*table_)[pos_];
	}

	virtual void mark_energy_computed() {
		(*computed_)[ pos_ ] = true;
	}

	virtual void mark_energy_uncomputed() {
		(*computed_)[ pos_ ] = false;
	}

	virtual bool energy_computed() const {
		return (*computed_)[ pos_ ];
	}

private:
	Size pos_;
	ScoreType score_type_;
	utility::vector1< Real > * table_;
	utility::vector1< bool > * computed_;
};


///////////////////////////////////////////////////////

class PeptideBondedNeighborConstIterator : public ResidueNeighborConstIterator
{
public:
	virtual ~PeptideBondedNeighborConstIterator(){}

	PeptideBondedNeighborConstIterator(
		Size const pos_in,
		ScoreType const st,
		utility::vector1< Real > const * table_in,
		utility::vector1< bool > const * computed_in
	):
		pos_( pos_in ),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & src ) {
		assert( dynamic_cast< PeptideBondedNeighborConstIterator const * >( &src ) );
		PeptideBondedNeighborConstIterator const & my_src( static_cast< PeptideBondedNeighborConstIterator const & >( src ) );
		pos_ = my_src.pos_;
		table_ = my_src.table_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborConstIterator const & operator ++ () {
		++pos_;
		return *this;
	}

	virtual bool operator == ( ResidueNeighborConstIterator const & other ) const {
		return ( residue_iterated_on() == other.residue_iterated_on() );
	}

	virtual bool operator != ( ResidueNeighborConstIterator const & other ) const {
		return !( *this == other );
	}

	virtual Size upper_neighbor_id() const {
		return pos_+1;
	}

	virtual Size lower_neighbor_id() const {
		return pos_;
	}

	virtual Size residue_iterated_on() const {
		return pos_;
	}

	virtual Size neighbor_id() const {
		return pos_+1;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const {
		emap[ score_type_ ] = (*table_)[pos_];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const {
		emap[ score_type_ ] += (*table_)[pos_];
	}

	virtual bool energy_computed() const {
		return (*computed_)[ pos_ ];
	}

private:
	Size pos_;
	ScoreType score_type_;
	utility::vector1< Real > const * table_;
	utility::vector1< bool > const * computed_;
};

///////////////////////////////////////////////////////////////////////////

class PeptideBondedEnergyContainer : public LREnergyContainer {
public:
	virtual
	~PeptideBondedEnergyContainer(){};

	virtual
	LREnergyContainerOP clone() const {
		return new PeptideBondedEnergyContainer( *this );
	}

	PeptideBondedEnergyContainer( Size const size_in, ScoreType const score_type_in ):
		size_( size_in ),
		score_type_( score_type_in ),
		table_( size_in, 0.0 ),
		computed_( size_in, false )
	{}

	virtual
	bool empty() const {
		return ( size_ == 0 );
	}

	virtual
	void
	set_num_nodes( Size size_in ) {
		size_ = size_in;
		table_.clear(); table_.resize( size_ , 0.0 );
		computed_.clear(); computed_.resize( size_,  false );
	}

	Size
	size() const {
		return size_;
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const {
		return new PeptideBondedNeighborConstIterator( resid-1, score_type_, &table_, &computed_ );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const {
		return new PeptideBondedNeighborConstIterator( std::min( resid+1, (int)size_ ), score_type_, &table_, &computed_ );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const {
		return new PeptideBondedNeighborConstIterator( resid, score_type_, &table_, &computed_ );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const {
		return const_neighbor_iterator_end( resid );
	}

	//////////////////// non-const versions
	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) {
		return new PeptideBondedNeighborIterator( resid-1, score_type_, &table_, &computed_ );
	}

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) {
		return new PeptideBondedNeighborIterator( std::min( resid+1, (int)size_ ), score_type_, &table_, &computed_ );
 	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid )
	{
		return new PeptideBondedNeighborIterator( resid, score_type_, &table_, &computed_ );
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) {
		return neighbor_iterator_end( resid );
	}

private:
	Size size_;
	ScoreType score_type_;

	utility::vector1< Real > table_;
	utility::vector1< bool > computed_;

};

} // namespace scoring
} // namespace core

#endif
