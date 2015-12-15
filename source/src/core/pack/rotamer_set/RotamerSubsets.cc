// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerSet/RotamerSets.cc
/// @brief  RotamerSets class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/rotamer_set/RotamerSubsets.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSubset.hh>

// Project headers
#include <core/conformation/Residue.hh>


// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
//#include <utility/integer_mapping.hh>
// C++


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

RotamerSubsets::RotamerSubsets(
	FixbbRotamerSets & source,
	utility::vector0< int > const & rotamer_subset
) :
	nmoltenres_( source.nmoltenres() ),
	total_residue_( source.total_residue() ),
	nrotamers_( rotamer_subset.size() ),
	set_of_rotamer_sets_( nmoltenres_ ),
	resid_2_moltenres_( source.resid_2_moltenres_vector() ),
	moltenres_2_resid_( source.moltenres_2_resid_vector() ),
	nrotamer_offsets_( nmoltenres_, 0 ),
	moltenres_for_rotamer_( nrotamers_, 0 ),
	nrotamers_for_moltenres_( nmoltenres_, 0 )
{
	for ( Size ii = 0; ii < rotamer_subset.size(); ++ii ) {
		Size const ii_moltres = source.moltenres_for_rotamer( rotamer_subset[ ii ] );
		++nrotamers_for_moltenres_[ ii_moltres ];
		moltenres_for_rotamer_[ ii + 1 ] = ii_moltres;
	}

	utility::vector1< utility::vector1< Size > > moltenres_subsets( nmoltenres_ );
	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		moltenres_subsets[ ii ].reserve( nrotamers_for_moltenres_[ ii ] );
		if ( ii > 1 ) {
			nrotamer_offsets_[ ii ] = nrotamer_offsets_[ ii - 1 ] + nrotamers_for_moltenres_[ ii - 1 ];
		}
	}

	for ( Size ii = 0; ii < rotamer_subset.size(); ++ii ) {
		Size const ii_rot = rotamer_subset[ ii ];
		Size const ii_moltres_id = source.moltenres_for_rotamer( ii_rot );
		Size const ii_local_rotno = source.rotid_on_moltenresidue( ii_rot );
		moltenres_subsets[ ii_moltres_id ].push_back( ii_local_rotno );
	}

	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		set_of_rotamer_sets_[ ii ] = RotamerSetOP( new RotamerSubset(
			*source.rotamer_set_for_moltenresidue( ii ),
			moltenres_subsets[ ii ]
			) );
	}
}


RotamerSubsets::~RotamerSubsets() {}

void
RotamerSubsets::update_offset_data()
{
	// count rotamers
	nrotamers_ = 0;
	for ( uint ii = 1; ii <= nmoltenres_; ++ii ) {
		nrotamers_for_moltenres_[ ii ] = set_of_rotamer_sets_[ii]->num_rotamers();
		nrotamers_ += set_of_rotamer_sets_[ii]->num_rotamers();
		if ( ii > 1 ) { nrotamer_offsets_[ ii ] = nrotamer_offsets_[ii - 1] + set_of_rotamer_sets_[ii - 1]->num_rotamers(); }
		else { nrotamer_offsets_[ ii ] = 0; }
	}

	moltenres_for_rotamer_.resize( nrotamers_ );
	uint count_rots_for_moltenres = 1;
	uint count_moltenres = 1;
	for ( uint ii = 1; ii <= nrotamers_; ++ii ) {
		moltenres_for_rotamer_[ ii ] = count_moltenres;
		if ( count_rots_for_moltenres == nrotamers_for_moltenres_[ count_moltenres ] ) {
			count_rots_for_moltenres = 1;
			++count_moltenres;
		} else {
			++count_rots_for_moltenres;
		}
	}
}


uint RotamerSubsets::nrotamers() const { return nrotamers_;}
uint RotamerSubsets::nrotamers_for_moltenres( uint mresid ) const
{
	return rotamer_set_for_moltenresidue( mresid )->num_rotamers();
}

uint RotamerSubsets::nmoltenres() const { return nmoltenres_;}

uint RotamerSubsets::total_residue() const { return total_residue_; }

uint
RotamerSubsets::moltenres_2_resid( uint mresid ) const { return moltenres_2_resid_[ mresid ]; }

uint
RotamerSubsets::resid_2_moltenres( uint resid ) const { return resid_2_moltenres_[ resid ]; }

uint
RotamerSubsets::moltenres_for_rotamer( uint rotid ) const { return moltenres_for_rotamer_[ rotid ]; }

uint
RotamerSubsets::res_for_rotamer( uint rotid ) const { return moltenres_2_resid( moltenres_for_rotamer( rotid ) ); }

core::conformation::ResidueCOP
RotamerSubsets::rotamer( uint rotid ) const
{
	return rotamer_set_for_residue( res_for_rotamer( rotid ) )->rotamer( rotid_on_moltenresidue( rotid ) );
}

core::conformation::ResidueCOP
RotamerSubsets::rotamer_for_moltenres( uint moltenres_id, uint rotamerid ) const
{
	return rotamer_set_for_residue( moltenres_id )->rotamer( rotamerid );
}


uint
RotamerSubsets::nrotamer_offset_for_moltenres( uint mresid ) const { return nrotamer_offsets_[ mresid ]; }

RotamerSetCOP
RotamerSubsets::rotamer_set_for_residue( uint resid ) const { return set_of_rotamer_sets_[ resid_2_moltenres( resid ) ]; }

RotamerSetOP
RotamerSubsets::rotamer_set_for_residue( uint resid )  { return set_of_rotamer_sets_[ resid_2_moltenres( resid ) ]; }

RotamerSetCOP
RotamerSubsets::rotamer_set_for_moltenresidue( uint moltenresid ) const { return set_of_rotamer_sets_[ moltenresid ]; }


RotamerSetOP
RotamerSubsets::rotamer_set_for_moltenresidue( uint moltenresid ) { return set_of_rotamer_sets_[ moltenresid ]; }

/// convert rotid in full rotamer enumeration into rotamer id on its source residue
uint
RotamerSubsets::rotid_on_moltenresidue( uint rotid ) const
{
	return rotid - nrotamer_offsets_[ moltenres_for_rotamer_[ rotid ] ];
}

/// convert moltenres rotid to id in full rotamer enumeration
uint
RotamerSubsets::moltenres_rotid_2_rotid( uint moltenres, uint moltenresrotid ) const
{
	return moltenresrotid + nrotamer_offsets_[ moltenres ];
}


void
RotamerSubsets::show( std::ostream & out ) const {
	out << "RotamerSubsets with " << nmoltenres_ << " molten residues for " << total_residue_ << " total residues and " << nrotamers_ << " rotamers." << std::endl;
	for ( core::Size ii(1); ii <= set_of_rotamer_sets_.size(); ++ii ) {
		out << ii << ": " << *(set_of_rotamer_sets_[ii]) << std::endl;
	}
}

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::pack::rotamer_set::RotamerSubsets::RotamerSubsets() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerSubsets::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::rotamer_set::FixbbRotamerSets >( this ) );
	arc( CEREAL_NVP( nmoltenres_ ) ); // uint
	arc( CEREAL_NVP( total_residue_ ) ); // uint
	arc( CEREAL_NVP( nrotamers_ ) ); // uint
	arc( CEREAL_NVP( set_of_rotamer_sets_ ) ); // RotamerSetVector
	arc( CEREAL_NVP( resid_2_moltenres_ ) ); // utility::vector1<uint>
	arc( CEREAL_NVP( moltenres_2_resid_ ) ); // utility::vector1<uint>
	arc( CEREAL_NVP( nrotamer_offsets_ ) ); // utility::vector1<uint>
	arc( CEREAL_NVP( moltenres_for_rotamer_ ) ); // utility::vector1<uint>
	arc( CEREAL_NVP( nrotamers_for_moltenres_ ) ); // utility::vector1<uint>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerSubsets::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::rotamer_set::FixbbRotamerSets >( this ) );
	arc( nmoltenres_ ); // uint
	arc( total_residue_ ); // uint
	arc( nrotamers_ ); // uint
	arc( set_of_rotamer_sets_ ); // RotamerSetVector
	arc( resid_2_moltenres_ ); // utility::vector1<uint>
	arc( moltenres_2_resid_ ); // utility::vector1<uint>
	arc( nrotamer_offsets_ ); // utility::vector1<uint>
	arc( moltenres_for_rotamer_ ); // utility::vector1<uint>
	arc( nrotamers_for_moltenres_ ); // utility::vector1<uint>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::RotamerSubsets );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::RotamerSubsets )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSubsets )
#endif // SERIALIZATION
