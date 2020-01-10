// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBPoseMap.cc
/// @brief  class to allow querying for pose resid with pdb chain/seqpos
/// @author Steven Lewis

// Unit headers
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// Project headers

// STL headers
#include <map>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.PDBPoseMap" );


/// @brief default constructor
PDBPoseMap::PDBPoseMap() :
	Super()
{}


/// @brief PDBInfo constructor
PDBPoseMap::PDBPoseMap( PDBInfo const & info ) :
	Super()
{
	fill( info );
}


/// @brief copy constructor
PDBPoseMap::PDBPoseMap( PDBPoseMap const & ) = default;


/// @brief default destructor
PDBPoseMap::~PDBPoseMap() = default;


/// @brief copy assignment
PDBPoseMap &
PDBPoseMap::operator =( PDBPoseMap const & m )
{
	if ( this != &m ) {
		pdb2pose_ = m.pdb2pose_;
		pdb2pose_noSegmentIDs_ = m.pdb2pose_noSegmentIDs_;
	}
	return *this;
}


Size
PDBPoseMap::find(
	char const chain,
	int const pdb_res,
	char const ins_code, /*= ' ', */
	std::string const & segmentID /* = "    " */
) const
{
	ResidueKey const rk( chain, pdb_res, ins_code, segmentID );
	auto const i = pdb2pose_.find( rk );

	if ( i == pdb2pose_.end() ) { // not found
		auto const j = pdb2pose_noSegmentIDs_.find( rk );
		if ( j == pdb2pose_noSegmentIDs_.end() ) { // not found again
			return 0;
		}
		return j->second;
	}
	return i->second; // return pose numbering
}


/// @brief insert pdb -> pose number mapping
/// @param[in] chain  chain id
/// @param[in] pdb_res  pdb residue numbering
/// @param[in] ins_code insertion code, use ' ' if no insertion code
/// @param[in] pose_res  pose numbering for residue
/// @remarks if the chain is equal to the PDBInfo's empty record character,
///  the insertion will be skipped
void
PDBPoseMap::insert(
	char const chain,
	int const pdb_res,
	char const ins_code,
	std::string const & segmentID,
	Size const pose_res
)
{
	if ( chain != PDBInfo::empty_record() ) {
		ResidueKey const rk( chain, pdb_res, ins_code, segmentID );
		pdb2pose_[ rk ] = pose_res;
		pdb2pose_noSegmentIDs_[ rk ] = pose_res;
	}
}



bool
PDBPoseMap::conditional_erase(
	char const chain,
	int const pdb_res,
	char const ins_code,
	std::string const & segmentID,
	Size const pose_res
) {
	ResidueKey const rk( chain, pdb_res, ins_code, segmentID );
	auto i = pdb2pose_.find( rk );
	if ( i != pdb2pose_.end() && i->second == pose_res ) {
		pdb2pose_.erase( i );
		pdb2pose_noSegmentIDs_.erase( rk );
		return true;
	}

	return false;
}


void
PDBPoseMap::erase(
	char const chain,
	int const pdb_res,
	char const ins_code,
	std::string const & segmentID /*= "    " */
)
{
	ResidueKey const rk( chain, pdb_res, ins_code, segmentID );
	pdb2pose_.erase( rk );
	pdb2pose_noSegmentIDs_.erase( rk );
}


/// @brief fill with corresponding pdb -> pose residue mapping
/// @note does not clear any currently existing mapping data
void
PDBPoseMap::fill( PDBInfo const & info )
{
	for ( Size i = 1; i <= info.nres(); ++i ) {
		insert( info.chain( i ), info.number( i ), info.icode( i ), info.segmentID( i ), i );
	}
}


} // pose
} // core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::PDBPoseMap::ResidueKey::save( Archive & arc ) const {
	arc( CEREAL_NVP( chainID ) ); // char
	arc( CEREAL_NVP( resSeq ) ); // int
	arc( CEREAL_NVP( iCode ) ); // char
	arc( CEREAL_NVP( segmentID ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::PDBPoseMap::ResidueKey::load( Archive & arc ) {
	arc( chainID ); // char
	arc( resSeq ); // int
	arc( iCode ); // char
	arc( segmentID ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::PDBPoseMap::ResidueKey );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::PDBPoseMap::save( Archive & arc ) const {
	arc( CEREAL_NVP( pdb2pose_ ) ); // Pdb2Pose
	arc( CEREAL_NVP( pdb2pose_noSegmentIDs_ ) ); // Pdb2Pose
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::PDBPoseMap::load( Archive & arc ) {
	arc( pdb2pose_ ); // Pdb2Pose
	arc( pdb2pose_noSegmentIDs_ ); // Pdb2Pose
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::PDBPoseMap );
CEREAL_REGISTER_TYPE( core::pose::PDBPoseMap )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_PDBPoseMap )
#endif // SERIALIZATION
