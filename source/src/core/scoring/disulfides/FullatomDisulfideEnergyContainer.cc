// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/ConstraintsEnergyContainer.cc
/// @brief  Constraints Energy Container class implementation
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/disulfides/FullatomDisulfideEnergyContainer.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/chemical/VariantType.hh>

// STL Headers
#include <utility/assert.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace disulfides {

DisulfResNeighbIterator::DisulfResNeighbIterator(
	FullatomDisulfideEnergyContainer * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

DisulfResNeighbIterator::DisulfResNeighbIterator(
	FullatomDisulfideEnergyContainer * owner
) :
	owner_( owner ),
	focused_residue_( FullatomDisulfideEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( FullatomDisulfideEnergyContainer::NO_DISULFIDE )
{}

DisulfResNeighbIterator::~DisulfResNeighbIterator() = default;

ResidueNeighborIterator &
DisulfResNeighbIterator::operator = ( ResidueNeighborIterator const & rhs)
{
	debug_assert( &(dynamic_cast< DisulfResNeighbIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;
}

// Incrementing an iterator in a list with exactly one element moves that
// iterator off the end of the list.
ResidueNeighborIterator const &
DisulfResNeighbIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	focused_residue_ = FullatomDisulfideEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = FullatomDisulfideEnergyContainer::NO_DISULFIDE;
	return *this;
}

bool
DisulfResNeighbIterator::operator == ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfResNeighbIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}

bool
DisulfResNeighbIterator::operator != ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfResNeighbIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}


Size
DisulfResNeighbIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

Size
DisulfResNeighbIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}

Size
DisulfResNeighbIterator::residue_iterated_on() const
{
	return focused_residue_;
}

Size
DisulfResNeighbIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}


void
DisulfResNeighbIterator::save_energy( EnergyMap const & emap )
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->save_energy( disulfide_index_, emap );
}


void
DisulfResNeighbIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}


void
DisulfResNeighbIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

void DisulfResNeighbIterator::mark_energy_computed()
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_computed( disulfide_index_ );
}

void DisulfResNeighbIterator::mark_energy_uncomputed()
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_uncomputed( disulfide_index_ );
}


bool
DisulfResNeighbIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}

////////////////////////////////////////////////////////////////////////
///// Disulfide Residue Neighbor Constant Iterator class implementation
////////////////////////////////////////////////////////////////////////

DisulfResNeighbConstIterator::DisulfResNeighbConstIterator(
	FullatomDisulfideEnergyContainer const * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

DisulfResNeighbConstIterator::DisulfResNeighbConstIterator(
	FullatomDisulfideEnergyContainer const * owner
) :
	owner_( owner ),
	focused_residue_( FullatomDisulfideEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( FullatomDisulfideEnergyContainer::NO_DISULFIDE )
{}

DisulfResNeighbConstIterator::~DisulfResNeighbConstIterator() = default;

ResidueNeighborConstIterator &
DisulfResNeighbConstIterator::operator = ( ResidueNeighborConstIterator const & rhs )
{
	debug_assert( &(dynamic_cast< DisulfResNeighbConstIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbConstIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;

}

ResidueNeighborConstIterator const &
DisulfResNeighbConstIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	focused_residue_ = FullatomDisulfideEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = FullatomDisulfideEnergyContainer::NO_DISULFIDE;
	return *this;
}

/// @brief returns true if the two edge-list iterators are equal
bool
DisulfResNeighbConstIterator::operator == ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfResNeighbConstIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbConstIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}


/// @brief returns true if the two edge-list iterators are not equal
bool
DisulfResNeighbConstIterator::operator != ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< DisulfResNeighbConstIterator const & > ( rhs )) );
	auto const & drni_rhs = static_cast< DisulfResNeighbConstIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}

Size
DisulfResNeighbConstIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

Size
DisulfResNeighbConstIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}


Size
DisulfResNeighbConstIterator::residue_iterated_on() const
{
	return focused_residue_;
}

Size
DisulfResNeighbConstIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}

/// @brief overwrites the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not zero out the other positions in the emap.
void
DisulfResNeighbConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}

/// @brief accumulates the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not touch the other positions in the emap.
void
DisulfResNeighbConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

bool
DisulfResNeighbConstIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != FullatomDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}


/////////////////////////////////////////////////////
/// Disulfide Energy Container Class Implementation
/////////////////////////////////////////////////////

Size const FullatomDisulfideEnergyContainer::NO_DISULFIDE( 0 );


FullatomDisulfideEnergyContainer::FullatomDisulfideEnergyContainer() = default;

bool
FullatomDisulfideEnergyContainer::empty() const
{
	return num_disulfides() == 0;
}


FullatomDisulfideEnergyContainer::FullatomDisulfideEnergyContainer( pose::Pose const & pose )
{
	find_disulfides( pose );
}

void
FullatomDisulfideEnergyContainer::update( pose::Pose const & pose )
{
	if ( disulfides_changed( pose ) ) find_disulfides( pose );
}

FullatomDisulfideEnergyContainer::~FullatomDisulfideEnergyContainer() = default;

LREnergyContainerOP
FullatomDisulfideEnergyContainer::clone() const
{
	FullatomDisulfideEnergyContainerOP dec( new FullatomDisulfideEnergyContainer );
	if ( !empty() ) {
		dec->disulfide_atom_indices_ = disulfide_atom_indices_;
		dec->disulfide_residue_types_ = disulfide_residue_types_;
		dec->resid_2_disulfide_index_ = resid_2_disulfide_index_;
		dec->disulfide_partners_ =  disulfide_partners_;
		dec->disulfide_info_ = disulfide_info_;
	}
	return dec;
}

void
FullatomDisulfideEnergyContainer::set_num_nodes( Size newsize )
{
	resid_2_disulfide_index_.resize( newsize, NO_DISULFIDE );
}

bool
FullatomDisulfideEnergyContainer::any_neighbors_for_residue( int resid ) const
{
	return resid <= (int) resid_2_disulfide_index_.size() && resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

bool
FullatomDisulfideEnergyContainer::any_upper_neighbors_for_residue( int resid ) const
{
	if ( resid <= (int) resid_2_disulfide_index_.size() ) {
		if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
			// disulfide_partners_ stores an ordered pair s.t. the first index is to the
			// lower residue and thes second index is to the upper residue; if we want to
			// know whether resid has an upper neighbor, what we need to know is if
			// it's the lower residue.
			return disulfide_partners_[ resid_2_disulfide_index_[ resid ] ].first == (Size) resid;
		}
	}
	return false;
}


ResidueNeighborConstIteratorOP
FullatomDisulfideEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );
	if ( resid <= (int) resid_2_disulfide_index_.size() && resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
FullatomDisulfideEnergyContainer::const_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this ) );
}

ResidueNeighborConstIteratorOP
FullatomDisulfideEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );
	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid <= (int)resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
FullatomDisulfideEnergyContainer::const_upper_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new DisulfResNeighbConstIterator( this ) );
}

ResidueNeighborIteratorOP
FullatomDisulfideEnergyContainer::neighbor_iterator_begin( int resid )
{
	debug_assert( !empty() );
	if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this ) );
	}
}

ResidueNeighborIteratorOP
FullatomDisulfideEnergyContainer::neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this ) );
}

ResidueNeighborIteratorOP
FullatomDisulfideEnergyContainer::upper_neighbor_iterator_begin( int resid )
{

	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid <= (int)resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this ) );
	}
}

ResidueNeighborIteratorOP
FullatomDisulfideEnergyContainer::upper_neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new DisulfResNeighbIterator( this ) );
}

bool
FullatomDisulfideEnergyContainer::disulfide_bonded( Size res1id, Size res2id ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ res1id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res2id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res1id ] == resid_2_disulfide_index_[ res2id ];
}

bool
FullatomDisulfideEnergyContainer::residue_forms_disulfide( Size resid ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

Size
FullatomDisulfideEnergyContainer::other_neighbor_id( Size resid ) const
{
	return other_neighbor_id( resid_2_disulfide_index_[ resid ], resid );
}


// Mutators
void
FullatomDisulfideEnergyContainer::save_energy( Size disulfide_index, EnergyMap const & emap )
{
	disulfide_info_[ disulfide_index ].first.dslf_ss_dst() = emap[ dslf_ss_dst ];
	disulfide_info_[ disulfide_index ].first.dslf_cs_ang() = emap[ dslf_cs_ang ];
	disulfide_info_[ disulfide_index ].first.dslf_ss_dih() = emap[ dslf_ss_dih ];
	disulfide_info_[ disulfide_index ].first.dslf_ca_dih() = emap[ dslf_ca_dih ];
	disulfide_info_[ disulfide_index ].first.dslf_cbs_ds() = emap[ dslf_cbs_ds ];
	disulfide_info_[ disulfide_index ].first.dslf_fa13()   = emap[ dslf_fa13 ];
}

void
FullatomDisulfideEnergyContainer::mark_energy_computed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = true;
}

void
FullatomDisulfideEnergyContainer::mark_energy_uncomputed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = false;
}

// Accessors
Size FullatomDisulfideEnergyContainer::lower_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].first;
}

Size FullatomDisulfideEnergyContainer::upper_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].second;
}

Size FullatomDisulfideEnergyContainer::other_neighbor_id( Size disulfide_index, Size resid ) const
{
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_partners_[ disulfide_index ].second :
		disulfide_partners_[ disulfide_index ].first );
}

DisulfideAtomIndices const &
FullatomDisulfideEnergyContainer::disulfide_atom_indices( Size resid ) const
{
	Size const disulfide_index( resid_2_disulfide_index_[ resid ] );
	debug_assert( disulfide_index != NO_DISULFIDE );
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_atom_indices_[ disulfide_index ].first :
		disulfide_atom_indices_[ disulfide_index ].second );
}


DisulfideAtomIndices const &
FullatomDisulfideEnergyContainer::other_neighbor_atom_indices( Size resid ) const
{
	Size const disulfide_index( resid_2_disulfide_index_[ resid ] );
	debug_assert( disulfide_index != NO_DISULFIDE );
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_atom_indices_[ disulfide_index ].second :
		disulfide_atom_indices_[ disulfide_index ].first );
}


void FullatomDisulfideEnergyContainer::accumulate_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslf_ss_dst ] += disulfide_info_[ disulfide_index ].first.dslf_ss_dst();
	emap[ dslf_cs_ang ] += disulfide_info_[ disulfide_index ].first.dslf_cs_ang();
	emap[ dslf_ss_dih ] += disulfide_info_[ disulfide_index ].first.dslf_ss_dih();
	emap[ dslf_ca_dih ] += disulfide_info_[ disulfide_index ].first.dslf_ca_dih();
	emap[ dslf_cbs_ds ] += disulfide_info_[ disulfide_index ].first.dslf_cbs_ds();
	emap[ dslf_fa13   ] += disulfide_info_[ disulfide_index ].first.dslf_fa13();
}

void FullatomDisulfideEnergyContainer::retrieve_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslf_ss_dst ] = disulfide_info_[ disulfide_index ].first.dslf_ss_dst();
	emap[ dslf_cs_ang ] = disulfide_info_[ disulfide_index ].first.dslf_cs_ang();
	emap[ dslf_ss_dih ] = disulfide_info_[ disulfide_index ].first.dslf_ss_dih();
	emap[ dslf_ca_dih ] = disulfide_info_[ disulfide_index ].first.dslf_ca_dih();
	emap[ dslf_cbs_ds ] = disulfide_info_[ disulfide_index ].first.dslf_cbs_ds();
	emap[ dslf_fa13   ] = disulfide_info_[ disulfide_index ].first.dslf_fa13();
}

bool FullatomDisulfideEnergyContainer::energy_computed( Size disulfide_index ) const
{
	return disulfide_info_[ disulfide_index ].second;
}


void
FullatomDisulfideEnergyContainer::find_disulfides( pose::Pose const & pose )
{
	Size nres = pose.size(), indep_res = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		indep_res = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		nres = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
	}

	disulfide_partners_.clear();
	disulfide_atom_indices_.clear();
	disulfide_info_.clear();
	resid_2_disulfide_index_.resize( nres );
	disulfide_residue_types_.resize( nres );
	std::fill( resid_2_disulfide_index_.begin(), resid_2_disulfide_index_.end(), NO_DISULFIDE );
	std::fill( disulfide_residue_types_.begin(), disulfide_residue_types_.end(), chemical::ResidueTypeCOP(nullptr) );


	Size count_disulfides( 0 );
	for ( Size ii = 1; ii <= indep_res; ++ii ) {
		if ( pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) &&
				resid_2_disulfide_index_[ ii ] == NO_DISULFIDE &&
				pose.residue_type( ii ).has( pose.residue_type( ii ).get_disulfide_atom_name() ) &&
				pose.residue_type( ii ).get_disulfide_atom_name() != "CEN" // full atom residue
				) {
			++count_disulfides;
			Size const ii_connect_atom = pose.residue( ii ).atom_index( pose.residue_type( ii ).get_disulfide_atom_name() );
			Size other_res( 0 );
			for ( Size jj = pose.residue( ii ).type().n_possible_residue_connections(); jj >= 1; --jj ) {
				if ( (Size) pose.residue( ii ).type().residue_connection( jj ).atomno() == ii_connect_atom ) {
					other_res = pose.residue( ii ).connect_map( jj ).resid();
					break;
				}
			}
			if ( other_res == 0 ) {
				std::cerr << "ERROR: Could not find disulfide partner for residue " << ii << std::endl;
				utility_exit();
			}
			debug_assert( other_res > ii );

			resid_2_disulfide_index_[ ii ] = count_disulfides;
			resid_2_disulfide_index_[ other_res ] = count_disulfides;
			disulfide_residue_types_[ ii ] = pose.residue_type_ptr( ii );
			disulfide_residue_types_[ other_res ] = pose.residue_type_ptr( other_res );
			disulfide_partners_.push_back( std::pair< Size, Size >( ii, other_res ) );
			disulfide_atom_indices_.push_back(
				std::pair< DisulfideAtomIndices, DisulfideAtomIndices > (
				DisulfideAtomIndices( pose.residue( ii ) ),
				DisulfideAtomIndices( pose.residue( other_res ) ) ));
			FullatomDisulfideEnergyComponents temp;
			disulfide_info_.push_back( std::pair< FullatomDisulfideEnergyComponents, bool > ( temp, false ) );
		}
	}

}

// we could do something like keep a flag for when minimization is occurring and assume that
// disulfide connectivity information does not change during the course of minimization...
bool
FullatomDisulfideEnergyContainer::disulfides_changed( pose::Pose const & pose )
{
	Size nres = pose.size(), indep_res = pose.size();
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		indep_res = core::pose::symmetry::symmetry_info( pose )->num_independent_residues();
		nres = core::pose::symmetry::symmetry_info( pose )->num_total_residues_without_pseudo();
	}
	Size const total_residue( indep_res );
	if ( resid_2_disulfide_index_.size() != nres ) return true;

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		if ( resid_2_disulfide_index_[ ii ] != NO_DISULFIDE ) {
			if ( disulfide_residue_types_[ ii ].get() != & ( pose.residue_type( ii ) )
					|| ! pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) /// subsumed by residue type check amw but reactivating
					|| pose.residue_type( ii ).get_disulfide_atom_name() == "CEN" //no longer full atom
					|| ! ( pose.residue( ii ).connect_map(
					pose.residue( ii ).type().residue_connection_id_for_atom(
					pose.residue( ii ).atom_index( pose.residue_type( ii ).get_disulfide_atom_name() ) ) ).resid()
					== other_neighbor_id( resid_2_disulfide_index_[ ii ], ii )
					) ) {
				return true;
			}
		} else if ( pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) ) {
			return true;
		}
	}
	return false;
}

Size FullatomDisulfideEnergyContainer::num_disulfides() const
{
	return disulfide_partners_.size();
}

Size FullatomDisulfideEnergyContainer::num_residues() const
{
	return resid_2_disulfide_index_.size();
}

}
}
}


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::disulfides::FullatomDisulfideEnergyContainer::save( Archive & arc ) const {

	arc( CEREAL_NVP( resid_2_disulfide_index_ ) ); // utility::vector1<Size>

	arc( CEREAL_NVP( disulfide_residue_types_ ) );

	arc( CEREAL_NVP( disulfide_partners_ ) ); // utility::vector1<std::pair<Size, Size> >
	arc( CEREAL_NVP( disulfide_atom_indices_ ) ); // utility::vector1<std::pair<DisulfideAtomIndices, DisulfideAtomIndices> >
	arc( CEREAL_NVP( disulfide_info_ ) ); // utility::vector1<std::pair<FullatomDisulfideEnergyComponents, _Bool> >
}

/// @details Deserialization method mostly auto-generated, but manually tweaked in order
/// to correctly deserialize ResidueTypeCOPs / resolve them to one of the globally-held
/// ResidueTypes.
template< class Archive >
void
core::scoring::disulfides::FullatomDisulfideEnergyContainer::load( Archive & arc ) {
	arc( resid_2_disulfide_index_ ); // utility::vector1<Size>

	arc( disulfide_residue_types_ );

	arc( disulfide_partners_ ); // utility::vector1<std::pair<Size, Size> >
	arc( disulfide_atom_indices_ ); // utility::vector1<std::pair<DisulfideAtomIndices, DisulfideAtomIndices> >
	arc( disulfide_info_ ); // utility::vector1<std::pair<FullatomDisulfideEnergyComponents, _Bool> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::disulfides::FullatomDisulfideEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::disulfides::FullatomDisulfideEnergyContainer )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::disulfides::FullatomDisulfideEnergyComponents::save( Archive & arc ) const {
	arc( CEREAL_NVP( dslf_ss_dst_ ) ); // Energy
	arc( CEREAL_NVP( dslf_cs_ang_ ) ); // Energy
	arc( CEREAL_NVP( dslf_ss_dih_ ) ); // Energy
	arc( CEREAL_NVP( dslf_ca_dih_ ) ); // Energy
	arc( CEREAL_NVP( dslf_cbs_ds_ ) ); // Energy
	arc( CEREAL_NVP( dslf_fa13_ ) ); // Energy
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::disulfides::FullatomDisulfideEnergyComponents::load( Archive & arc ) {
	arc( dslf_ss_dst_ ); // Energy
	arc( dslf_cs_ang_ ); // Energy
	arc( dslf_ss_dih_ ); // Energy
	arc( dslf_ca_dih_ ); // Energy
	arc( dslf_cbs_ds_ ); // Energy
	arc( dslf_fa13_ ); // Energy
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::disulfides::FullatomDisulfideEnergyComponents );

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_disulfides_FullatomDisulfideEnergyContainer )
#endif // SERIALIZATION
