// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/ConstraintsEnergyContainer.cc
/// @brief  Constraints Energy Container class implementation
/// @author Spencer Bliven <blivens@u.washington.edu>

// Unit headers
#include <core/scoring/disulfides/CentroidDisulfideEnergyContainer.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// STL Headers
#include <utility/assert.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.disulfides.CentroidDisulfideEnergyContainer" );

/// @brief constructor
CentroidDisulfideNeighborIterator::CentroidDisulfideNeighborIterator(
	CentroidDisulfideEnergyContainer * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

/// @brief constructor, default to no disulfide bond
CentroidDisulfideNeighborIterator::CentroidDisulfideNeighborIterator(
	CentroidDisulfideEnergyContainer * owner
)
:
	owner_( owner ),
	focused_residue_( CentroidDisulfideEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( CentroidDisulfideEnergyContainer::NO_DISULFIDE )
{}

CentroidDisulfideNeighborIterator::~CentroidDisulfideNeighborIterator()
{}

/// @brief Assignment
ResidueNeighborIterator const &
CentroidDisulfideNeighborIterator::operator = ( ResidueNeighborIterator const & rhs)
{
	debug_assert( &(dynamic_cast< CentroidDisulfideNeighborIterator const & > ( rhs )) );
	CentroidDisulfideNeighborIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;
}

/// @note Incrementing an iterator in a list with exactly one element moves that
/// iterator off the end of the list.
ResidueNeighborIterator const &
CentroidDisulfideNeighborIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	focused_residue_ = CentroidDisulfideEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = CentroidDisulfideEnergyContainer::NO_DISULFIDE;
	return *this;
}

bool
CentroidDisulfideNeighborIterator::operator == ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< CentroidDisulfideNeighborIterator const & > ( rhs )) );
	CentroidDisulfideNeighborIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}

bool
CentroidDisulfideNeighborIterator::operator != ( ResidueNeighborIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< CentroidDisulfideNeighborIterator const & > ( rhs )) );
	CentroidDisulfideNeighborIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}

/// @brief Get the higher-numbered residue for this disulfide bond
Size
CentroidDisulfideNeighborIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

/// @brief Get the lower-numbered residue for this disulfide bond
Size
CentroidDisulfideNeighborIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}

/// @brief Which residue are we looking for disulfide bonds to?
Size
CentroidDisulfideNeighborIterator::residue_iterated_on() const
{
	return focused_residue_;
}

/// @brief Which residue is disulfide bonded to the current residue?
Size
CentroidDisulfideNeighborIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}

/// @brief Save the specified energies for this disulfide to the
/// CentroidDisulfideEnergyContainer associated with this iterator.
void
CentroidDisulfideNeighborIterator::save_energy( EnergyMap const & emap )
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->save_energy( disulfide_index_, emap );
}

/// @brief Get the energies for the current disulfide bond from the
/// CentroidDisulfideEnergyContainer associated with this iterator.
void
CentroidDisulfideNeighborIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}

/// @brief Add some energies to the totals already in CentroidDisulfideEnergyContainer
void
CentroidDisulfideNeighborIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

void CentroidDisulfideNeighborIterator::mark_energy_computed()
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_computed( disulfide_index_ );
}

void CentroidDisulfideNeighborIterator::mark_energy_uncomputed()
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->mark_energy_uncomputed( disulfide_index_ );
}


bool
CentroidDisulfideNeighborIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}

////////////////////////////////////////////////////////////////////////
///// Disulfide Residue Neighbor Constant Iterator class implementation
////////////////////////////////////////////////////////////////////////

CentroidDisulfideNeighborConstIterator::CentroidDisulfideNeighborConstIterator(
	CentroidDisulfideEnergyContainer const * owner,
	Size focused_residue,
	Size disulfide_index
) :
	owner_( owner ),
	focused_residue_( focused_residue ),
	disulfide_index_( disulfide_index )
{}

CentroidDisulfideNeighborConstIterator::CentroidDisulfideNeighborConstIterator(
	CentroidDisulfideEnergyContainer const * owner
) :
	owner_( owner ),
	focused_residue_( CentroidDisulfideEnergyContainer::NO_DISULFIDE ),
	disulfide_index_( CentroidDisulfideEnergyContainer::NO_DISULFIDE )
{}

CentroidDisulfideNeighborConstIterator::~CentroidDisulfideNeighborConstIterator()
{}

ResidueNeighborConstIterator const &
CentroidDisulfideNeighborConstIterator::operator = ( ResidueNeighborConstIterator const & rhs )
{
	debug_assert( &(dynamic_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs )) );
	CentroidDisulfideNeighborConstIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs );

	owner_ = drni_rhs.owner_;
	focused_residue_ = drni_rhs.focused_residue_;
	disulfide_index_ = drni_rhs.disulfide_index_;
	return *this;

}

ResidueNeighborConstIterator const &
CentroidDisulfideNeighborConstIterator::operator ++ ()
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	focused_residue_ = CentroidDisulfideEnergyContainer::NO_DISULFIDE;
	disulfide_index_ = CentroidDisulfideEnergyContainer::NO_DISULFIDE;
	return *this;
}

/// @brief returns true if the two edge-list iterators are equal
bool
CentroidDisulfideNeighborConstIterator::operator == ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs )) );
	CentroidDisulfideNeighborConstIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs );

	return ( owner_ == drni_rhs.owner_ &&
		focused_residue_ == drni_rhs.focused_residue_ &&
		disulfide_index_ == drni_rhs.disulfide_index_ );
}


/// @brief returns true if the two edge-list iterators are not equal
bool
CentroidDisulfideNeighborConstIterator::operator != ( ResidueNeighborConstIterator const & rhs ) const
{
	debug_assert( &( dynamic_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs )) );
	CentroidDisulfideNeighborConstIterator const & drni_rhs = static_cast< CentroidDisulfideNeighborConstIterator const & > ( rhs );
	return ( owner_ != drni_rhs.owner_ ||
		focused_residue_ != drni_rhs.focused_residue_ ||
		disulfide_index_ != drni_rhs.disulfide_index_ );
}

Size
CentroidDisulfideNeighborConstIterator::upper_neighbor_id() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->upper_neighbor_id( disulfide_index_ );
}

Size
CentroidDisulfideNeighborConstIterator::lower_neighbor_id() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->lower_neighbor_id( disulfide_index_ );
}


Size
CentroidDisulfideNeighborConstIterator::residue_iterated_on() const
{
	return focused_residue_;
}

Size
CentroidDisulfideNeighborConstIterator::neighbor_id() const
{
	return owner_->other_neighbor_id( disulfide_index_, focused_residue_ );
}

/// @brief overwrites the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not zero out the other positions in the emap.
void
CentroidDisulfideNeighborConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->retrieve_energy( disulfide_index_, emap );
}

/// @brief accumulates the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not touch the other positions in the emap.
void
CentroidDisulfideNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	owner_->accumulate_energy( disulfide_index_, emap );
}

bool
CentroidDisulfideNeighborConstIterator::energy_computed() const
{
	debug_assert( disulfide_index_ != CentroidDisulfideEnergyContainer::NO_DISULFIDE );
	return owner_->energy_computed( disulfide_index_ );
}


/////////////////////////////////////////////////////
/// Disulfide Energy Container Class Implementation
/////////////////////////////////////////////////////

Size const CentroidDisulfideEnergyContainer::NO_DISULFIDE( 0 );


CentroidDisulfideEnergyContainer::CentroidDisulfideEnergyContainer()
{}

bool
CentroidDisulfideEnergyContainer::empty() const
{
	return num_disulfides() == 0;
}


CentroidDisulfideEnergyContainer::CentroidDisulfideEnergyContainer( pose::Pose const & pose )
{
	find_disulfides( pose );
}

void
CentroidDisulfideEnergyContainer::update( pose::Pose const & pose )
{
	if ( disulfides_changed( pose ) ) find_disulfides( pose );
}

CentroidDisulfideEnergyContainer::~CentroidDisulfideEnergyContainer()
{}

LREnergyContainerOP
CentroidDisulfideEnergyContainer::clone() const
{
	CentroidDisulfideEnergyContainerOP dec( new CentroidDisulfideEnergyContainer );
	if ( !empty() ) {
		dec->disulfide_atom_indices_ = disulfide_atom_indices_;
		dec->disulfide_residue_types_ = disulfide_residue_types_;
		dec->resid_2_disulfide_index_ = resid_2_disulfide_index_;
		dec->disulfide_partners_ =  disulfide_partners_;
		dec->disulfide_info_ = disulfide_info_;
	}
	return dec;
}

bool
CentroidDisulfideEnergyContainer::any_neighbors_for_residue( int resid ) const
{
	return (Size) resid <= resid_2_disulfide_index_.size() && resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

bool
CentroidDisulfideEnergyContainer::any_upper_neighbors_for_residue( int resid ) const
{
	if ( (Size) resid <= resid_2_disulfide_index_.size() ) {
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
CentroidDisulfideEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );
	if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
CentroidDisulfideEnergyContainer::const_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this ) );
}

ResidueNeighborConstIteratorOP
CentroidDisulfideEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	debug_assert( !empty() );

	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid >= 0 &&
			Size( resid ) <= resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this ) );
	}
}

ResidueNeighborConstIteratorOP
CentroidDisulfideEnergyContainer::const_upper_neighbor_iterator_end( int ) const
{
	debug_assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CentroidDisulfideNeighborConstIterator( this ) );
}

ResidueNeighborIteratorOP
CentroidDisulfideEnergyContainer::neighbor_iterator_begin( int resid )
{
	debug_assert( !empty() );
	if ( resid_2_disulfide_index_[ resid ] != NO_DISULFIDE ) {
		return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this ) );
	}
}

ResidueNeighborIteratorOP
CentroidDisulfideEnergyContainer::neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this ) );
}

ResidueNeighborIteratorOP
CentroidDisulfideEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	// cppcheck flags this but it is fine -- the limits check is happening in the right order!
	if ( resid >= 0 &&
			Size( resid ) <= resid_2_disulfide_index_.size() &&
			resid_2_disulfide_index_[ resid ] != NO_DISULFIDE &&
			(Size) resid < other_neighbor_id( resid_2_disulfide_index_[ resid ], resid ) ) {
		return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this, resid, resid_2_disulfide_index_[ resid ]  ) );
	} else {
		return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this ) );
	}
}

ResidueNeighborIteratorOP
CentroidDisulfideEnergyContainer::upper_neighbor_iterator_end( int )
{
	debug_assert( !empty() );
	return ResidueNeighborIteratorOP( new CentroidDisulfideNeighborIterator( this ) );
}

bool
CentroidDisulfideEnergyContainer::disulfide_bonded( Size res1id, Size res2id ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ res1id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res2id ] != NO_DISULFIDE &&
		resid_2_disulfide_index_[ res1id ] == resid_2_disulfide_index_[ res2id ];
}

bool
CentroidDisulfideEnergyContainer::residue_forms_disulfide( Size resid ) const
{
	if ( empty() ) return false;
	return resid_2_disulfide_index_[ resid ] != NO_DISULFIDE;
}

Size
CentroidDisulfideEnergyContainer::other_neighbor_id( Size resid ) const
{
	return other_neighbor_id( resid_2_disulfide_index_[ resid ], resid );
}


// Mutators
void
CentroidDisulfideEnergyContainer::save_energy( Size disulfide_index, EnergyMap const & emap )
{
	disulfide_info_[ disulfide_index ].first.dslfc_cen_dst() = emap[ dslfc_cen_dst ];
	disulfide_info_[ disulfide_index ].first.dslfc_cb_dst() = emap[ dslfc_cb_dst ];
	disulfide_info_[ disulfide_index ].first.dslfc_ang() = emap[ dslfc_ang ];
	disulfide_info_[ disulfide_index ].first.dslfc_cb_dih() = emap[ dslfc_cb_dih ];
	disulfide_info_[ disulfide_index ].first.dslfc_bb_dih() = emap[ dslfc_bb_dih ];
}

void
CentroidDisulfideEnergyContainer::mark_energy_computed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = true;
}

void
CentroidDisulfideEnergyContainer::mark_energy_uncomputed( Size disulfide_index )
{
	disulfide_info_[ disulfide_index ].second = false;
}

// Accessors
Size CentroidDisulfideEnergyContainer::lower_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].first;
}

Size CentroidDisulfideEnergyContainer::upper_neighbor_id( Size disulfide_index ) const
{
	return disulfide_partners_[ disulfide_index ].second;
}

Size CentroidDisulfideEnergyContainer::other_neighbor_id( Size disulfide_index, Size resid ) const
{
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_partners_[ disulfide_index ].second :
		disulfide_partners_[ disulfide_index ].first );
}

DisulfideAtomIndices const &
CentroidDisulfideEnergyContainer::disulfide_atom_indices( Size resid ) const
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
CentroidDisulfideEnergyContainer::other_neighbor_atom_indices( Size resid ) const
{
	Size const disulfide_index( resid_2_disulfide_index_[ resid ] );
	debug_assert( disulfide_index != NO_DISULFIDE );
	debug_assert( disulfide_partners_[ disulfide_index ].first == resid ||
		disulfide_partners_[ disulfide_index ].second == resid );
	return ( resid == disulfide_partners_[ disulfide_index ].first ?
		disulfide_atom_indices_[ disulfide_index ].second :
		disulfide_atom_indices_[ disulfide_index ].first );
}


void CentroidDisulfideEnergyContainer::accumulate_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslfc_cen_dst ] += disulfide_info_[ disulfide_index ].first.dslfc_cen_dst();
	emap[ dslfc_cb_dst  ] += disulfide_info_[ disulfide_index ].first.dslfc_cb_dst();
	emap[ dslfc_ang     ] += disulfide_info_[ disulfide_index ].first.dslfc_ang();
	emap[ dslfc_cb_dih  ] += disulfide_info_[ disulfide_index ].first.dslfc_cb_dih();
	emap[ dslfc_bb_dih  ] += disulfide_info_[ disulfide_index ].first.dslfc_bb_dih();
}

void CentroidDisulfideEnergyContainer::retrieve_energy( Size disulfide_index, EnergyMap & emap ) const
{
	emap[ dslfc_cen_dst ] = disulfide_info_[ disulfide_index ].first.dslfc_cen_dst();
	emap[ dslfc_cb_dst  ] = disulfide_info_[ disulfide_index ].first.dslfc_cb_dst();
	emap[ dslfc_ang     ] = disulfide_info_[ disulfide_index ].first.dslfc_ang();
	emap[ dslfc_cb_dih  ] = disulfide_info_[ disulfide_index ].first.dslfc_cb_dih();
	emap[ dslfc_bb_dih  ] = disulfide_info_[ disulfide_index ].first.dslfc_bb_dih();
}

bool CentroidDisulfideEnergyContainer::energy_computed( Size disulfide_index ) const
{
	return disulfide_info_[ disulfide_index ].second;
}


void
CentroidDisulfideEnergyContainer::find_disulfides( pose::Pose const & pose )
{
	TR.Debug << "In find_disulfides():" << std::endl;

	disulfide_partners_.clear();
	disulfide_atom_indices_.clear();
	disulfide_info_.clear();
	resid_2_disulfide_index_.resize( pose.total_residue() );
	disulfide_residue_types_.resize( pose.total_residue() );
	std::fill( resid_2_disulfide_index_.begin(), resid_2_disulfide_index_.end(), NO_DISULFIDE );
	std::fill( disulfide_residue_types_.begin(), disulfide_residue_types_.end(), chemical::ResidueTypeCOP(0) );

	Size count_disulfides( 0 );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue res = pose.residue( ii );
		if ( res.type().is_disulfide_bonded() &&
				//res.aa() == chemical::aa_cys &&
				res.has_variant_type( chemical::DISULFIDE ) &&
				resid_2_disulfide_index_[ ii ] == NO_DISULFIDE &&
				pose.residue_type( ii ).has( "CEN" )
				) {
			++count_disulfides;
			Size const ii_connect_atom( res.atom_index( "CEN" ) );
			Size other_res_ii( 0 );
			for ( Size jj = 1; jj <= res.type().n_residue_connections(); ++jj ) {
				if ( (Size) res.type().residue_connection( jj ).atomno() == ii_connect_atom ) {
					other_res_ii = res.connect_map( jj ).resid();
					break;
				}
			}
			if ( other_res_ii == 0 ) {
				TR.Error << "ERROR: Could not find disulfide partner for residue " << ii << std::endl;
				utility_exit();
			}
			debug_assert( other_res_ii > ii );
			//Can only bond residues of the same residue type set (eg centroid to centroid)
			debug_assert( pose.residue_type(other_res_ii).residue_type_set().name() ==
				pose.residue_type(ii).residue_type_set().name() );

			TR.Debug << "Found disulf between " << ii << " and " << other_res_ii << std::endl;

			resid_2_disulfide_index_[ ii ] = count_disulfides;
			resid_2_disulfide_index_[ other_res_ii ] = count_disulfides;
			disulfide_residue_types_[ ii ] = pose.residue_type( ii ).get_self_ptr();
			disulfide_residue_types_[ other_res_ii ] = pose.residue_type( other_res_ii ).get_self_ptr();
			disulfide_partners_.push_back( std::pair< Size, Size >( ii, other_res_ii ) );
			disulfide_atom_indices_.push_back( std::pair< DisulfideAtomIndices, DisulfideAtomIndices > (
				DisulfideAtomIndices( pose.residue(ii ) ), DisulfideAtomIndices( pose.residue( other_res_ii ) ) ));
			CentroidDisulfideEnergyComponents temp;
			disulfide_info_.push_back( std::pair< CentroidDisulfideEnergyComponents, bool > ( temp, false ) );

			debug_assert(! empty());
		}
	}
	TR.Debug << "Found " << num_disulfides() << " DS" << std::endl;
}

// we could do something like keep a flag for when minimization is occurring and assume that
// disulfide connectivity information does not change during the course of minimization...
bool
CentroidDisulfideEnergyContainer::disulfides_changed( pose::Pose const & pose )
{
	Size const total_residue( pose.total_residue() );
	if ( resid_2_disulfide_index_.size() != total_residue ) return true;

	for ( Size ii = 1; ii <= total_residue; ++ii ) {
		if ( resid_2_disulfide_index_[ ii ] != NO_DISULFIDE ) {
			conformation::Residue res = pose.residue( ii );
			if ( //res.aa() != chemical::aa_cys ||
					!res.type().is_disulfide_bonded() ||
					disulfide_residue_types_[ ii ].get() != & (pose.residue_type( ii )) ||
					/// subsumed by residue type check ! pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) ||
					! pose.residue_type( ii ).has( "CEN" ) || // not centroid
					res.connect_map(
					res.type().residue_connection_id_for_atom(
					res.atom_index( "CEN" ) ) ).resid() !=
					other_neighbor_id( resid_2_disulfide_index_[ ii ], ii ) ) {
				return true;
			}
		} else if ( pose.residue( ii ).aa() == chemical::aa_cys &&
				pose.residue( ii ).has_variant_type( chemical::DISULFIDE ) ) {
			return true;
		}
	}
	return false;
}

Size CentroidDisulfideEnergyContainer::num_disulfides() const
{
	return disulfide_partners_.size();
}


}
}
}

