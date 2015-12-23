// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/DenseEnergyContainer.cc
/// @brief  A container for storing all-against-all energies as the upper triangle of a matrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/DenseEnergyContainer.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray2D.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

DenseNeighborIterator::~DenseNeighborIterator(){}

DenseNeighborIterator::DenseNeighborIterator(
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

ResidueNeighborIterator &
DenseNeighborIterator::operator = ( ResidueNeighborIterator const & src )
{
	debug_assert( dynamic_cast< DenseNeighborIterator const * >( &src ) );
	DenseNeighborIterator const & my_src( static_cast< DenseNeighborIterator const & >( src ) );
	pos1_ = my_src.pos1_;
	pos2_ = my_src.pos2_;
	table_ = my_src.table_;
	computed_ = my_src.computed_;
	return *this;
}

ResidueNeighborIterator const &
DenseNeighborIterator::operator ++ ()
{
	++pos2_;
	if ( pos2_ == pos1_ ) ++pos2_;
	return *this;
}

bool
DenseNeighborIterator::operator == ( ResidueNeighborIterator const & other ) const
{
	return ( residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool
DenseNeighborIterator::operator != ( ResidueNeighborIterator const & other ) const
{
	return !( *this == other );
}

Size
DenseNeighborIterator::upper_neighbor_id() const
{
	return std::max( pos1_, pos2_ );
}

Size
DenseNeighborIterator::lower_neighbor_id() const
{
	return std::min( pos1_, pos2_ );
}

Size
DenseNeighborIterator::residue_iterated_on() const
{
	return pos1_;
}

Size
DenseNeighborIterator::neighbor_id() const
{
	return pos2_;
}

void
DenseNeighborIterator::save_energy( EnergyMap const & emap )
{
	Real const energy( emap[ score_type_ ] );
	(*table_)( pos1_, pos2_ ) = energy;
	(*table_)( pos2_, pos1_ ) = energy;
}

void
DenseNeighborIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] = (*table_)(pos1_,pos2_);
}

void
DenseNeighborIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] += (*table_)(pos1_, pos2_);
}

void
DenseNeighborIterator::mark_energy_computed()
{
	(*computed_)( pos1_, pos2_ ) = true;
	(*computed_)( pos2_, pos1_ ) = true;
}

void
DenseNeighborIterator::mark_energy_uncomputed()
{
	(*computed_)( pos1_, pos2_ ) = false;
	(*computed_)( pos2_, pos1_ ) = false;
}

bool
DenseNeighborIterator::energy_computed() const
{
	return (*computed_)( pos1_, pos2_ );
}

/////////////////////////////////////////////////////

DenseNeighborConstIterator::~DenseNeighborConstIterator(){}

DenseNeighborConstIterator::DenseNeighborConstIterator(
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

ResidueNeighborConstIterator &
DenseNeighborConstIterator::operator = ( ResidueNeighborConstIterator const & src )
{
	debug_assert( dynamic_cast< DenseNeighborConstIterator const * >( &src ) );
	DenseNeighborConstIterator const & my_src( static_cast< DenseNeighborConstIterator const & >( src ) );
	pos1_ = my_src.pos1_;
	pos2_ = my_src.pos2_;
	table_ = my_src.table_;
	computed_ = my_src.computed_;
	return *this;
}

ResidueNeighborConstIterator const &
DenseNeighborConstIterator::operator ++ ()
{
	++pos2_;
	if ( pos2_ == pos1_ ) ++pos2_;
	return *this;
}

bool
DenseNeighborConstIterator::operator == ( ResidueNeighborConstIterator const & other ) const
{
	return ( residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool
DenseNeighborConstIterator::operator != ( ResidueNeighborConstIterator const & other ) const
{
	return !( *this == other );
}

Size
DenseNeighborConstIterator::upper_neighbor_id() const
{
	return std::max( pos1_, pos2_ );
}

Size
DenseNeighborConstIterator::lower_neighbor_id() const
{
	return std::min( pos1_, pos2_ );
}

Size
DenseNeighborConstIterator::residue_iterated_on() const
{
	return pos1_;
}

Size
DenseNeighborConstIterator::neighbor_id() const
{
	return pos2_;
}

void
DenseNeighborConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] = (*table_)(pos1_,pos2_);
}

void
DenseNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ score_type_ ] += (*table_)(pos1_, pos2_);
}

bool
DenseNeighborConstIterator::energy_computed() const
{
	return (*computed_)( pos1_, pos2_ );
}

/////////////////////////////////////////////////////////////////////////

DenseEnergyContainer::~DenseEnergyContainer() {}

LREnergyContainerOP
DenseEnergyContainer::clone() const
{
	return LREnergyContainerOP( new DenseEnergyContainer( *this ) );
}

DenseEnergyContainer::DenseEnergyContainer( Size const size_in, ScoreType const score_type_in ):
	size_( size_in ),
	score_type_( score_type_in ),
	table_( size_, size_, 0.0 ),
	computed_( size_, size_, false )
{}

bool
DenseEnergyContainer::empty() const
{
	return ( size_ == 0 );
}

void
DenseEnergyContainer::set_num_nodes( Size size_in ) {
	size_ = size_in;
	table_.dimension( size_ , size_, 0.0 );
	computed_.dimension( size_, size_, false );
}

bool
DenseEnergyContainer::any_neighbors_for_residue( int ) const {
	return true;
}

bool
DenseEnergyContainer::any_upper_neighbors_for_residue( int resid ) const
{
	return (Size) resid < size_;
}

Size
DenseEnergyContainer::size() const
{
	return size_;
}

//////////////////// const versions

ResidueNeighborConstIteratorOP
DenseEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, resid==1 ? 2 : 1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborConstIteratorOP
DenseEnergyContainer::const_neighbor_iterator_end( int resid ) const
{
	return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, size_ + 1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborConstIteratorOP
DenseEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	return ResidueNeighborConstIteratorOP( new DenseNeighborConstIterator( resid, resid+1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborConstIteratorOP
DenseEnergyContainer::const_upper_neighbor_iterator_end( int resid ) const
{
	return const_neighbor_iterator_end( resid );
}

//////////////////// non-const versions

ResidueNeighborIteratorOP
DenseEnergyContainer::neighbor_iterator_begin( int resid )
{
	return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, 1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborIteratorOP
DenseEnergyContainer::neighbor_iterator_end( int resid )
{
	return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, size_ + 1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborIteratorOP
DenseEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	return ResidueNeighborIteratorOP( new DenseNeighborIterator( resid, resid+1, score_type_, &table_, &computed_ ) );
}

ResidueNeighborIteratorOP
DenseEnergyContainer::upper_neighbor_iterator_end( int resid )
{
	return neighbor_iterator_end( resid );
}

} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::DenseEnergyContainer::DenseEnergyContainer() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::DenseEnergyContainer::save( Archive & arc ) const {
	arc( CEREAL_NVP( size_ ) ); // Size
	arc( CEREAL_NVP( score_type_ ) ); // enum core::scoring::ScoreType
	arc( CEREAL_NVP( table_ ) ); // ObjexxFCL::FArray2D<Real>
	arc( CEREAL_NVP( computed_ ) ); // ObjexxFCL::FArray2D<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::DenseEnergyContainer::load( Archive & arc ) {
	arc( size_ ); // Size
	arc( score_type_ ); // enum core::scoring::ScoreType
	arc( table_ ); // ObjexxFCL::FArray2D<Real>
	arc( computed_ ); // ObjexxFCL::FArray2D<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::DenseEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::DenseEnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_DenseEnergyContainer )
#endif // SERIALIZATION
