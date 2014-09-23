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
	virtual ~PeptideBondedNeighborIterator(){}

	PeptideBondedNeighborIterator(
		Size const base_in,
		Size const pos_in,
		utility::vector1< ScoreType > const st,
		utility::vector1< utility::vector1< Real > > * table_in,
		utility::vector1< bool > * computed_in
	):
		base_( base_in ),
		pos_( pos_in ),
		score_types_( st ),
		tables_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & src ) {
		assert( dynamic_cast< PeptideBondedNeighborIterator const * >( &src ) );
		PeptideBondedNeighborIterator const & my_src( static_cast< PeptideBondedNeighborIterator const & >( src ) );
		base_ = my_src.base_;
		pos_ = my_src.pos_;
		tables_ = my_src.tables_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborIterator const & operator ++ () {
		++pos_;
		if (pos_ == base_) ++pos_;
		return *this;
	}

	virtual bool operator == ( ResidueNeighborIterator const & other ) const
	{
		return ( residue_iterated_on() == other.residue_iterated_on() &&
						 neighbor_id() == other.neighbor_id() );
	}

	virtual bool operator != ( ResidueNeighborIterator const & other ) const
	{
		return !( *this == other );
	}

	virtual Size upper_neighbor_id() const {
		return std::max(pos_,base_);
	}

	virtual Size lower_neighbor_id() const {
		return std::min(pos_,base_);
	}

	virtual Size residue_iterated_on() const {
		return base_;
	}

	virtual Size neighbor_id() const {
		return pos_;
	}

	virtual void save_energy( EnergyMap const & emap ) {
		for (Size i=1; i<=score_types_.size(); ++i)
		{
			Real const energy( emap[ score_types_[i] ] );
			(*tables_)[ std::min(pos_,base_) ][i] = energy;
		}
	}

	virtual void retrieve_energy( EnergyMap & emap ) const {
		for (Size i=1; i<=score_types_.size(); ++i)
			emap[ score_types_[i] ] = (*tables_)[std::min(pos_,base_)][i];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const {
		for (Size i=1; i<=score_types_.size(); ++i)
			emap[ score_types_[i] ] += (*tables_)[std::min(pos_,base_)][i];
	}

	virtual void mark_energy_computed() {
		(*computed_)[ std::min(pos_,base_) ] = true;
	}

	virtual void mark_energy_uncomputed() {
		(*computed_)[ std::min(pos_,base_) ] = false;
	}

	virtual bool energy_computed() const {
		return (*computed_)[ std::min(pos_,base_) ];
	}

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
	virtual ~PeptideBondedNeighborConstIterator(){}

	PeptideBondedNeighborConstIterator(
		Size const base_in,
		Size const pos_in,
		utility::vector1< ScoreType > const st,
		utility::vector1< utility::vector1< Real > > const * table_in,
		utility::vector1< bool > const * computed_in
	):
		base_( base_in ),
		pos_( pos_in ),
		score_types_( st ),
		tables_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & src ) {
		assert( dynamic_cast< PeptideBondedNeighborConstIterator const * >( &src ) );
		PeptideBondedNeighborConstIterator const & my_src( static_cast< PeptideBondedNeighborConstIterator const & >( src ) );
		pos_ = my_src.pos_;
		tables_ = my_src.tables_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborConstIterator const & operator ++ () {
		++pos_;
		if (pos_ == base_) ++pos_;
		return *this;
	}

	virtual bool operator == ( ResidueNeighborConstIterator const & other ) const {
		return ( residue_iterated_on() == other.residue_iterated_on() &&
						 neighbor_id() == other.neighbor_id() );
	}

	virtual bool operator != ( ResidueNeighborConstIterator const & other ) const {
		return !( *this == other );
	}

	virtual Size upper_neighbor_id() const {
		return std::max(pos_,base_);
	}

	virtual Size lower_neighbor_id() const {
		return std::min(pos_,base_);
	}

	virtual Size residue_iterated_on() const {
		return base_;
	}

	virtual Size neighbor_id() const {
		return pos_;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const {
		for (Size i=1; i<=score_types_.size(); ++i)
			emap[ score_types_[i] ] = (*tables_)[std::min(pos_,base_)][i];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const {
		for (Size i=1; i<=score_types_.size(); ++i)
			emap[ score_types_[i] ] += (*tables_)[std::min(pos_,base_)][i];
	}

	virtual bool energy_computed() const {
		return (*computed_)[ std::min(pos_,base_) ];
	}

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
	~PeptideBondedEnergyContainer(){};

	virtual
	LREnergyContainerOP clone() const {
		return LREnergyContainerOP( new PeptideBondedEnergyContainer( *this ) );
	}

	PeptideBondedEnergyContainer( Size const size_in, utility::vector1< ScoreType > const score_type_in, Size offset_in=0 ):
		size_( size_in ),
		offset_( offset_in ),
		score_types_( score_type_in ),
		computed_( size_in+offset_in, false )
	{
		int nscoretypes = score_type_in.size();
		tables_.resize( size_in+offset_in, utility::vector1< core::Real >(nscoretypes,0) );
	}

	virtual
	bool empty() const {
		return ( size_ == 0 );
	}

	using LREnergyContainer::set_num_nodes;

	virtual
	void
	set_num_nodes( Size size_in, Size offset_in=0 ) {
		size_ = size_in;
		offset_ = offset_in;
		int nscoretypes = score_types_.size();
		tables_.clear();
		tables_.resize( size_+offset_, utility::vector1< core::Real >(nscoretypes,0) );
		computed_.clear();
		computed_.resize( size_+offset_,  false );
	}

	Size
	size() const {
		return size_;
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const {
		int beginat = std::min( resid-1, (int)(offset_+size_+1) );
		if (resid==(int)(offset_+1)) beginat = offset_+2;
		if (resid<=(int)(offset_)) beginat = 1; // sometimes arises in symmetry
		if (resid> (int)(offset_+size_)) beginat = 1; // sometimes arises in symmetry
		return ResidueNeighborConstIteratorOP( new PeptideBondedNeighborConstIterator( resid, beginat, score_types_, &tables_, &computed_ ) );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const {
		int endat = std::min( resid+2, (int)(offset_+size_+1) );
		if (resid<=(int)(offset_)) endat = 1; // sometimes arises in symmetry
		if (resid>(int)(offset_+size_)) endat = 1; // sometimes arises in symmetry
		return ResidueNeighborConstIteratorOP( new PeptideBondedNeighborConstIterator( resid, endat, score_types_, &tables_, &computed_ ) );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const {
		int beginat = std::min( resid+1, (int)(offset_+size_+1) );
		if (resid<=(int)(offset_)) beginat = 1; // sometimes arises in symmetry
		if (resid> (int)(offset_+size_)) beginat = 1; // sometimes arises in symmetry
		return ResidueNeighborConstIteratorOP( new PeptideBondedNeighborConstIterator( resid, beginat, score_types_, &tables_, &computed_ ) );
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
		int beginat = std::min( resid-1, (int)(offset_+size_+1) );
		if (resid==(int)(offset_+1)) beginat = offset_+2;
		if (resid<=(int)(offset_)) beginat = 1; // sometimes arises in symmetry
		if (resid> (int)(offset_+size_)) beginat = 1; // sometimes arises in symmetry
		return ResidueNeighborIteratorOP( new PeptideBondedNeighborIterator( resid, beginat, score_types_, &tables_, &computed_ ) );
	}

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) {
		int endat = std::min( resid+2, (int)(offset_+size_+1) );
		if (resid<=(int)(offset_)) endat = 1; // sometimes arises in symmetry
		if (resid>(int)(offset_+size_)) endat = 1; // sometimes arises in symmetry
		return ResidueNeighborIteratorOP( new PeptideBondedNeighborIterator( resid, endat, score_types_, &tables_, &computed_ ) );
 	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid )
	{
		int beginat = std::min( resid+1, (int)(offset_+size_+1) );
		if (resid<=(int)(offset_)) beginat = 1; // sometimes arises in symmetry
		if (resid> (int)(offset_+size_)) beginat = 1; // sometimes arises in symmetry
		return ResidueNeighborIteratorOP( new PeptideBondedNeighborIterator( resid, beginat, score_types_, &tables_, &computed_ ) );
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) {
		return neighbor_iterator_end( resid );
	}

private:
	Size size_, offset_;
	utility::vector1< ScoreType > score_types_;
	utility::vector1< utility::vector1< Real > > tables_;
	utility::vector1< bool > computed_;

};

} // namespace scoring
} // namespace core

#endif
