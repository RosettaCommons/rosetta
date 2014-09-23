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


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class OneToAllNeighborIterator : public ResidueNeighborIterator
{
public:
	virtual ~OneToAllNeighborIterator(){}

	OneToAllNeighborIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > * table_in,
		utility::vector1< bool > * computed_in
	):
		pos1_( pos1_in ),
		pos2_( pos2_in ),
		operating_on_pos1_(operating_on_pos1_in),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & src )
	{
		assert( dynamic_cast< OneToAllNeighborIterator const * >( &src ) );
		OneToAllNeighborIterator const & my_src( static_cast< OneToAllNeighborIterator const & >( src ) );
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
		return pos1_;  // "upper" is always the fixed res, which is pos1_
	}

	virtual Size lower_neighbor_id() const
	{
		return pos2_;
	}

	virtual Size residue_iterated_on() const
	{
		return operating_on_pos1_? pos1_ : pos2_;
	}

	virtual Size neighbor_id() const
	{
		return operating_on_pos1_? pos2_ : pos1_;
	}

	virtual void save_energy( EnergyMap const & emap )
	{
		Real const energy( emap[ score_type_ ] );
		(*table_)[ pos2_ ] = energy;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] = (*table_)[pos2_];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] += (*table_)[pos2_];
	}

	virtual void mark_energy_computed()
	{
		//std::cerr << "OO  mark_energy_computed( " << pos2_ << " )\n";
		(*computed_)[ pos2_ ] = true;
	}

	virtual void mark_energy_uncomputed()
	{
		//std::cerr << "XX  mark_energy_uncomputed( " << pos2_ << " )\n";
		(*computed_)[ pos2_ ] = false;
	}

	virtual bool energy_computed() const
	{
		return (*computed_)[ pos2_ ];
	}

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
	virtual ~OneToAllNeighborConstIterator(){}

	OneToAllNeighborConstIterator(
		Size const pos1_in,
		Size const pos2_in,
		bool const operating_on_pos1_in,
		ScoreType const st,
		utility::vector1< Real > const * table_in,
		utility::vector1< bool > const * computed_in
	):
		pos1_( pos1_in ),
		pos2_( pos2_in ),
		operating_on_pos1_(operating_on_pos1_in),
		score_type_( st ),
		table_( table_in ),
		computed_( computed_in )
	{}

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & src )
	{
		assert( dynamic_cast< OneToAllNeighborConstIterator const * >( &src ) );
		OneToAllNeighborConstIterator const & my_src( static_cast< OneToAllNeighborConstIterator const & >( src ) );
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
		return pos1_;  // "upper" is always the fixed res, which is pos1_
	}

	virtual Size lower_neighbor_id() const
	{
		return pos2_;
	}

	virtual Size residue_iterated_on() const
	{
		return operating_on_pos1_? pos1_ : pos2_;
	}

	virtual Size neighbor_id() const
	{
		return operating_on_pos1_? pos2_ : pos1_;
	}

	virtual void retrieve_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] = (*table_)[pos2_];
	}

	virtual void accumulate_energy( EnergyMap & emap ) const
	{
		emap[ score_type_ ] += (*table_)[pos2_];
	}

	virtual bool energy_computed() const
	{
		return (*computed_)[ pos2_ ];
	}

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
	~OneToAllEnergyContainer(){};

	virtual
	LREnergyContainerOP clone() const
	{
		return LREnergyContainerOP( new OneToAllEnergyContainer( *this ) );
	}

	OneToAllEnergyContainer( int const fixed_res_idx, Size const size_in, ScoreType const score_type_in ):
		fixed_( fixed_res_idx ),
		size_( size_in ),
		score_type_( score_type_in ),
		table_( size_, 0.0 ),
		computed_( size_, false )
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
		table_.clear(); table_.resize( size_ , 0.0 );
		computed_.clear(); computed_.resize( size_,  false );
	}

	Size
	size() const
	{
		return size_;
	}

	int
	fixed() const
	{
		return fixed_;
	}

	//////////////////// const versions
	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const
	{
		if (resid == fixed_) {
			// loop over ALL tgts
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			//std::cerr << "START fixed " << fixed_ << " , resid " << resid << std::endl;
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , resid, false, score_type_, &table_, &computed_ ) );
		}
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const
	{
		if (resid == fixed_) {
			// loop over ALL tgts
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , size_ + 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			if (resid+1 == fixed_) {
				//std::cerr << "END fixed " << fixed_ << " , resid " << resid+2 << std::endl;
				return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , resid + 2, false, score_type_, &table_, &computed_ ) );
			} else  {
				//std::cerr << "END fixed " << fixed_ << " , resid " << resid+1 << std::endl;
				return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , resid + 1, false, score_type_, &table_, &computed_ ) );
			}
		}
	}

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const
	{
		if (resid == fixed_) {
			// loop over NOTHING
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , size_ + 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_ , resid, false, score_type_, &table_, &computed_ ) );
		}
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
		if (resid == fixed_) {
			// loop over ALL tgts
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , resid, false, score_type_, &table_, &computed_ ) );
		}
	}

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid )
	{
		if (resid == fixed_) {
			// loop over ALL tgts
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , size_ + 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			if (resid+1 == fixed_)
				return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , resid + 2, false, score_type_, &table_, &computed_ ) );
			else
				return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , resid + 1, false, score_type_, &table_, &computed_ ) );
		}
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid )
	{
		if (resid == fixed_) {
			// loop over NOTHING
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , size_ + 1, true, score_type_, &table_, &computed_ ) );
		} else {
			// loop over fixed only
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_ , resid, false, score_type_, &table_, &computed_ ) );
		}
	}

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid )
	{
		return neighbor_iterator_end( resid );
	}

private:
	int fixed_;
	Size /*const*/ size_;
	ScoreType score_type_;

	utility::vector1< Real > table_;
	utility::vector1< bool > computed_;

};

} // namespace scoring
} // namespace core

#endif
