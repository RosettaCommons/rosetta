// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/LongRangeEnergyContainer.hh
/// @brief  A container interface for storing and scoring long range energies
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



namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class DenseNeighborIterator : public ResidueNeighborIterator
{
public:
	virtual ~DenseNeighborIterator(){}

	DenseNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > * table_in,
		ObjexxFCL::FArray2D< bool > * computed_in
	):
		pos1_( pos1_in ),
		pos2_( pos2_in ),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & src )
	{
		assert( dynamic_cast< DenseNeighborIterator const * >( &src ) );
		DenseNeighborIterator const & my_src( static_cast< DenseNeighborIterator const & >( src ) );
		pos1_ = my_src.pos1_;
		pos2_ = my_src.pos2_;
		table_ = my_src.table_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborIterator const & operator ++ ()
	{
		++pos2_;
		if ( pos2_ == pos1_ ) ++pos2_;
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

	virtual Size upper_neighbor_id() const
	{
		return std::max( pos1_, pos2_ );
	}

	virtual Size lower_neighbor_id() const
	{
		return std::min( pos1_, pos2_ );
	}

	virtual Size residue_iterated_on() const
	{
		return pos1_;
	}

	virtual Size neighbor_id() const
	{
		return pos2_;
	}

	virtual void save_energy( EnergyMap const & emap )
	{
		Real const energy( emap[ score_type_ ] );
		(*table_)( pos1_, pos2_ ) = energy;
		(*table_)( pos2_, pos1_ ) = energy;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] = (*table_)(pos1_,pos2_);
	}

	virtual void accumulate_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] += (*table_)(pos1_, pos2_);
	}

	virtual void mark_energy_computed()
	{
		(*computed_)( pos1_, pos2_ ) = true;
		(*computed_)( pos2_, pos1_ ) = true;
	}

	virtual void mark_energy_uncomputed()
	{
		(*computed_)( pos1_, pos2_ ) = false;
		(*computed_)( pos2_, pos1_ ) = false;
	}

	virtual bool energy_computed() const
	{
		return (*computed_)( pos1_, pos2_ );
	}

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
public:
	virtual ~DenseNeighborConstIterator(){}

	DenseNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		ScoreType const st,
		ObjexxFCL::FArray2D< Real > const * table_in,
		ObjexxFCL::FArray2D< bool > const * computed_in
	):
		pos1_( pos1_in ),
		pos2_( pos2_in ),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & src )
	{
		assert( dynamic_cast< DenseNeighborConstIterator const * >( &src ) );
		DenseNeighborConstIterator const & my_src( static_cast< DenseNeighborConstIterator const & >( src ) );
		pos1_ = my_src.pos1_;
		pos2_ = my_src.pos2_;
		table_ = my_src.table_;
		computed_ = my_src.computed_;
		return *this;
	}

	virtual ResidueNeighborConstIterator const & operator ++ ()
	{
		++pos2_;
		if ( pos2_ == pos1_ ) ++pos2_;
		return *this;
	}

	virtual bool operator == ( ResidueNeighborConstIterator const & other ) const
	{
		return ( residue_iterated_on() == other.residue_iterated_on() &&
						 neighbor_id() == other.neighbor_id() );
	}

	virtual bool operator != ( ResidueNeighborConstIterator const & other ) const
	{
		return !( *this == other );
	}

	virtual Size upper_neighbor_id() const
	{
		return std::max( pos1_, pos2_ );
	}

	virtual Size lower_neighbor_id() const
	{
		return std::min( pos1_, pos2_ );
	}

	virtual Size residue_iterated_on() const
	{
		return pos1_;
	}

	virtual Size neighbor_id() const
	{
		return pos2_;
	}

// 	virtual void save_energy( EnergyMap const & emap )
// 	{
// 		Real const energy( emap[ score_type_ ] );
// 		(*table_)( pos1_, pos2_ ) = energy;
// 		(*table_)( pos2_, pos1_ ) = energy;
// 	}

	virtual void retrieve_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] = (*table_)(pos1_,pos2_);
	}

	virtual void accumulate_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] += (*table_)(pos1_, pos2_);
	}

// 	virtual void mark_energy_computed()
// 	{
// 		(*computed_)( pos1_, pos2_ ) = true;
// 		(*computed_)( pos2_, pos1_ ) = true;
// 	}

// 	virtual void mark_energy_uncomputed()
// 	{
// 		(*computed_)( pos1_, pos2_ ) = false;
// 		(*computed_)( pos2_, pos1_ ) = false;
// 	}

	virtual bool energy_computed() const
	{
		return (*computed_)( pos1_, pos2_ );
	}

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
	virtual
	~DenseEnergyContainer(){};

	virtual
	LREnergyContainerOP clone() const
	{
		return LREnergyContainerOP( new DenseEnergyContainer( *this ) );
	}

	DenseEnergyContainer( Size const size_in, ScoreType const score_type_in ):
		size_( size_in ),
		score_type_( score_type_in ),
		table_( size_, size_, 0.0 ),
		computed_( size_, size_, false )
	{}

	virtual
	bool empty() const
	{
		return ( size_ == 0 );
	}

	virtual
	void
	set_num_nodes( Size size_in ) {
		size_ = size_in;
		table_.dimension( size_ , size_, 0.0 );
		computed_.dimension( size_, size_, false );
	}

	Size
	size() const
	{
		return size_;
	}

	//////////////////// const versions
	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const
	{
		return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, resid==1 ? 2 : 1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const
	{
		return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, size_ + 1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const
	{
		return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, resid+1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const
	{
		return const_neighbor_iterator_end( resid );
	}

	//////////////////// non-const versions
	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid )
	{
		return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, 1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid )
	{
		return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, size_ + 1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid )
	{
		return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, resid+1, score_type_, &table_, &computed_ ) );
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid )
	{
		return neighbor_iterator_end( resid );
	}

private:
	Size /*const*/ size_;
	ScoreType score_type_;

	ObjexxFCL::FArray2D< Real > table_;
	ObjexxFCL::FArray2D< bool > computed_;

};

} // namespace scoring
} // namespace core

#endif
