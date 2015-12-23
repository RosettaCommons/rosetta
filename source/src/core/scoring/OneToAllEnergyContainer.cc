// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/OneToAllEnergyContainer.cc
/// @brief  A container interface for storing and scoring long range energies
/// @author Frank DiMaio

// Unit headers
#include <core/scoring/OneToAllEnergyContainer.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

/////////////////////////////////////////////////////

OneToAllNeighborIterator::~OneToAllNeighborIterator(){}

OneToAllNeighborIterator::OneToAllNeighborIterator(
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

ResidueNeighborIterator & OneToAllNeighborIterator::operator = ( ResidueNeighborIterator const & src )
{
	debug_assert( dynamic_cast< OneToAllNeighborIterator const * >( &src ) );
	OneToAllNeighborIterator const & my_src( static_cast< OneToAllNeighborIterator const & >( src ) );
	pos1_ = my_src.pos1_;
	pos2_ = my_src.pos2_;
	table_ = my_src.table_;
	computed_ = my_src.computed_;
	return *this;
}

ResidueNeighborIterator const & OneToAllNeighborIterator::operator ++ ()
{
	++pos2_;
	if ( pos2_ == pos1_ ) ++pos2_;
	return *this;
}

bool OneToAllNeighborIterator::operator == ( ResidueNeighborIterator const & other ) const
{
	return ( residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool OneToAllNeighborIterator::operator != ( ResidueNeighborIterator const & other ) const
{
	return !( *this == other );
}

Size OneToAllNeighborIterator::upper_neighbor_id() const
{
	return pos1_;  // "upper" is always the fixed res, which is pos1_
}

Size OneToAllNeighborIterator::lower_neighbor_id() const
{
	return pos2_;
}

Size OneToAllNeighborIterator::residue_iterated_on() const
{
	return operating_on_pos1_? pos1_ : pos2_;
}

Size OneToAllNeighborIterator::neighbor_id() const
{
	return operating_on_pos1_? pos2_ : pos1_;
}

void OneToAllNeighborIterator::save_energy( EnergyMap const & emap )
{
	Real const energy( emap[ score_type_ ] );
	(*table_)[ pos2_ ] = energy;
}

void OneToAllNeighborIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] = (*table_)[pos2_];
}

void OneToAllNeighborIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] += (*table_)[pos2_];
}

void OneToAllNeighborIterator::mark_energy_computed()
{
	//std::cerr << "OO  mark_energy_computed( " << pos2_ << " )\n";
	(*computed_)[ pos2_ ] = true;
}

void OneToAllNeighborIterator::mark_energy_uncomputed()
{
	//std::cerr << "XX  mark_energy_uncomputed( " << pos2_ << " )\n";
	(*computed_)[ pos2_ ] = false;
}

bool OneToAllNeighborIterator::energy_computed() const
{
	return (*computed_)[ pos2_ ];
}


/////////////////////////////////////////////////////

OneToAllNeighborConstIterator::~OneToAllNeighborConstIterator(){}

OneToAllNeighborConstIterator::OneToAllNeighborConstIterator(
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

ResidueNeighborConstIterator & OneToAllNeighborConstIterator::operator = ( ResidueNeighborConstIterator const & src )
{
	debug_assert( dynamic_cast< OneToAllNeighborConstIterator const * >( &src ) );
	OneToAllNeighborConstIterator const & my_src( static_cast< OneToAllNeighborConstIterator const & >( src ) );
	pos1_ = my_src.pos1_;
	pos2_ = my_src.pos2_;
	table_ = my_src.table_;
	computed_ = my_src.computed_;
	return *this;
}

ResidueNeighborConstIterator const & OneToAllNeighborConstIterator::operator ++ ()
{
	++pos2_;
	if ( pos2_ == pos1_ ) ++pos2_;
	return *this;
}

bool OneToAllNeighborConstIterator::operator == ( ResidueNeighborConstIterator const & other ) const
{
	return ( residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool OneToAllNeighborConstIterator::operator != ( ResidueNeighborConstIterator const & other ) const
{
	return !( *this == other );
}

Size OneToAllNeighborConstIterator::upper_neighbor_id() const
{
	return pos1_;  // "upper" is always the fixed res, which is pos1_
}

Size OneToAllNeighborConstIterator::lower_neighbor_id() const
{
	return pos2_;
}

Size OneToAllNeighborConstIterator::residue_iterated_on() const
{
	return operating_on_pos1_? pos1_ : pos2_;
}

Size OneToAllNeighborConstIterator::neighbor_id() const
{
	return operating_on_pos1_? pos2_ : pos1_;
}

void OneToAllNeighborConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] = (*table_)[pos2_];
}

void OneToAllNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] += (*table_)[pos2_];
}

bool OneToAllNeighborConstIterator::energy_computed() const
{
	return (*computed_)[ pos2_ ];
}


/////////////////////////////////////////////////////////////////////////

OneToAllEnergyContainer::~OneToAllEnergyContainer() {}

LREnergyContainerOP OneToAllEnergyContainer::clone() const
{
	return LREnergyContainerOP( new OneToAllEnergyContainer( *this ) );
}

OneToAllEnergyContainer::OneToAllEnergyContainer( int const fixed_res_idx, Size const size_in, ScoreType const score_type_in ):
	fixed_( fixed_res_idx ),
	size_( size_in ),
	score_type_( score_type_in ),
	table_( size_, 0.0 ),
	computed_( size_, false )
{}

bool OneToAllEnergyContainer::empty() const
{
	return ( size_ == 0 );
}

void
OneToAllEnergyContainer::set_num_nodes( Size size_in ) {
	size_ = size_in;
	table_.clear(); table_.resize( size_ , 0.0 );
	computed_.clear(); computed_.resize( size_,  false );
}

bool
OneToAllEnergyContainer::any_neighbors_for_residue( int ) const {
	return true;
}

bool
OneToAllEnergyContainer::any_upper_neighbors_for_residue( int resid ) const {
	return resid == fixed_;
}

Size
OneToAllEnergyContainer::size() const
{
	return size_;
}

int
OneToAllEnergyContainer::fixed() const
{
	return fixed_;
}

//////////////////// const versions
ResidueNeighborConstIteratorOP
OneToAllEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	if ( resid == fixed_ ) {
		// loop over ALL tgts
		return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		//std::cerr << "START fixed " << fixed_ << " , resid " << resid << std::endl;
		return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, resid, false, score_type_, &table_, &computed_ ) );
	}
}

ResidueNeighborConstIteratorOP
OneToAllEnergyContainer::const_neighbor_iterator_end( int resid ) const
{
	if ( resid == fixed_ ) {
		// loop over ALL tgts
		return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, size_ + 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		if ( resid+1 == fixed_ ) {
			//std::cerr << "END fixed " << fixed_ << " , resid " << resid+2 << std::endl;
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, resid + 2, false, score_type_, &table_, &computed_ ) );
		} else  {
			//std::cerr << "END fixed " << fixed_ << " , resid " << resid+1 << std::endl;
			return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, resid + 1, false, score_type_, &table_, &computed_ ) );
		}
	}
}

ResidueNeighborConstIteratorOP
OneToAllEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	if ( resid == fixed_ ) {
		// loop over NOTHING
		return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, size_ + 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		return ResidueNeighborConstIteratorOP( new OneToAllNeighborConstIterator( fixed_, resid, false, score_type_, &table_, &computed_ ) );
	}
}

ResidueNeighborConstIteratorOP
OneToAllEnergyContainer::const_upper_neighbor_iterator_end( int resid ) const
{
	return const_neighbor_iterator_end( resid );
}

//////////////////// non-const versions
ResidueNeighborIteratorOP
OneToAllEnergyContainer::neighbor_iterator_begin( int resid )
{
	if ( resid == fixed_ ) {
		// loop over ALL tgts
		return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, resid, false, score_type_, &table_, &computed_ ) );
	}
}

ResidueNeighborIteratorOP
OneToAllEnergyContainer::neighbor_iterator_end( int resid )
{
	if ( resid == fixed_ ) {
		// loop over ALL tgts
		return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, size_ + 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		if ( resid+1 == fixed_ ) {
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, resid + 2, false, score_type_, &table_, &computed_ ) );
		} else {
			return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, resid + 1, false, score_type_, &table_, &computed_ ) );
		}
	}
}

ResidueNeighborIteratorOP
OneToAllEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	if ( resid == fixed_ ) {
		// loop over NOTHING
		return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, size_ + 1, true, score_type_, &table_, &computed_ ) );
	} else {
		// loop over fixed only
		return ResidueNeighborIteratorOP( new OneToAllNeighborIterator( fixed_, resid, false, score_type_, &table_, &computed_ ) );
	}
}

ResidueNeighborIteratorOP
OneToAllEnergyContainer::upper_neighbor_iterator_end( int resid )
{
	return neighbor_iterator_end( resid );
}

} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::OneToAllEnergyContainer::OneToAllEnergyContainer() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::OneToAllEnergyContainer::save( Archive & arc ) const {
	arc( CEREAL_NVP( fixed_ ) ); // int
	arc( CEREAL_NVP( size_ ) ); // Size
	arc( CEREAL_NVP( score_type_ ) ); // enum core::scoring::ScoreType
	arc( CEREAL_NVP( table_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( computed_ ) ); // utility::vector1<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::OneToAllEnergyContainer::load( Archive & arc ) {
	arc( fixed_ ); // int
	arc( size_ ); // Size
	arc( score_type_ ); // enum core::scoring::ScoreType
	arc( table_ ); // utility::vector1<Real>
	arc( computed_ ); // utility::vector1<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::OneToAllEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::OneToAllEnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_OneToAllEnergyContainer )
#endif // SERIALIZATION
