// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/PDBPoseMap.cc
/// @brief  class to allow querying for pose resid with pdb chain/seqpos
/// @author Steven Lewis

// Unit headers
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// Project headers

// STL headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace pose {

static thread_local basic::Tracer TR( "core.pose.PDBPoseMap" );


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
PDBPoseMap::PDBPoseMap( PDBPoseMap const & map ) :
	Super( map ),
	pdb2pose_( map.pdb2pose_ )
{}


/// @brief default destructor
PDBPoseMap::~PDBPoseMap()
{}


/// @brief copy assignment
PDBPoseMap &
PDBPoseMap::operator =( PDBPoseMap const & m )
{
	if ( this != &m ) {
		pdb2pose_ = m.pdb2pose_;
	}
	return *this;
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
	Size const pose_res
)
{
	if ( chain != PDBInfo::empty_record() ) {
		pdb2pose_[ ResidueKey( chain, pdb_res, ins_code ) ] = pose_res;
	}
}


/// @brief fill with corresponding pdb -> pose residue mapping
/// @note does not clear any currently existing mapping data
void
PDBPoseMap::fill( PDBInfo const & info )
{
	for ( Size i = 1; i <= info.nres(); ++i ) {
		insert( info.chain( i ), info.number( i ), info.icode( i ), i );
	}
}


} // pose
} // core
