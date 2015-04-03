// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/EntityCorrespondence.cc
/// @brief  Implementation for class EntityCorrespondence
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/pack_daemon/EntityCorrespondence.hh>

// C/C++ headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>

namespace protocols {
namespace pack_daemon {

enum ec_funcs
{
	func_add_resid_to_entity_list = 1,
	func_entity_for_residue,
	func_n_residues_for_entity,
	func_residues_for_entity_begin,
	func_residues_for_entity_end
};

EntityCorrespondence::EntityCorrespondence() :
	funcnames_( func_residues_for_entity_end )
{
	funcnames_[ func_add_resid_to_entity_list  ] = "add_resid_to_entity_list";
	funcnames_[ func_entity_for_residue        ] = "entity_for_residue";
	funcnames_[ func_n_residues_for_entity     ] = "n_residues_for_entity";
	funcnames_[ func_residues_for_entity_begin ] = "residues_for_entity_begin";
	funcnames_[ func_residues_for_entity_end   ] = "residues_for_entity_end";

}
EntityCorrespondence::~EntityCorrespondence() {}

EntityCorrespondence::EntityCorrespondence(EntityCorrespondence const & src ) :
	parent(),
	pose_( src.pose_ ),
	pdb_pose_map_( src.pdb_pose_map_ ),
	entity_id_2_resids_( src.entity_id_2_resids_ ),
	resid_2_entity_id_( src.resid_2_entity_id_ ),
	funcnames_( src.funcnames_ )
{}

EntityCorrespondence const &
EntityCorrespondence::operator = ( EntityCorrespondence const & other )
{
	if ( this != &other ) {
		pose_ = other.pose_;
		pdb_pose_map_ = other.pdb_pose_map_;
		entity_id_2_resids_ = other.entity_id_2_resids_;
		resid_2_entity_id_ = other.resid_2_entity_id_;
		funcnames_ = other.funcnames_;
	}
	return *this;
}


void EntityCorrespondence::set_pose( core::pose::PoseCOP pose )
{
	using namespace core::pose;

	if ( pose_ == 0 || pose_->total_residue() != pose->total_residue() ) {
		resid_2_entity_id_.resize( pose->total_residue() );
		std::fill( resid_2_entity_id_.begin(), resid_2_entity_id_.end(), 0 );
	}
	pose_ = pose;
	if ( ! pose_->pdb_info() ) {
		PDBInfo info(*pose_);
		pdb_pose_map_ = core::pose::PDBPoseMapCOP( core::pose::PDBPoseMapOP( new PDBPoseMap( info ) ) ); //use PDBInfo ctor for PDBPoseMap
	} else {
		pdb_pose_map_ = core::pose::PDBPoseMapCOP( core::pose::PDBPoseMapOP( new PDBPoseMap( pose_->pdb_info()->pdb2pose()) ) );
	}
}

void EntityCorrespondence::set_num_entities( Size num_entities )
{
	entity_id_2_resids_.resize( num_entities );
	std::fill( resid_2_entity_id_.begin(), resid_2_entity_id_.end(), 0 );
}

/// @details File format for entity correspondence:
/// Column 1: entity id.  Columns 2 and 3 PDB id (resid+[optional insertion code] chain)
void
EntityCorrespondence::initialize_from_correspondence_file( std::istream & instream )
{
	Size count_correspondences( 0 );
	Size line_count( 0 );
	int PDBnum = 0;
	std::string PDBnum_string;
	char chain = '_';
	Size entity_id( 0 ), residue_index( 0 );
	while ( instream ) {
		char icode = ' ';
		++line_count;
		std::string line;
		getline( instream, line );
		if ( line.size() == 0 ) continue; // Ignore blank lines

		std::istringstream linestream( line );
		if (( linestream >> entity_id ).fail() ) {
			throw utility::excn::EXCN_Msg_Exception(
				"Failed to read entity id on line " + utility::to_string( line_count ) +
				" of the EntityCorrespondence file:\n" + line );
		}
		++count_correspondences;
		if ( entity_id > num_entities() ) {
			throw utility::excn::EXCN_Msg_Exception(
				"Entity ID read on line " + utility::to_string( line_count ) +
				" of the EntityCorrespondence file exceeds the number of entities: " +
				utility::to_string( entity_id ) + " vs " + utility::to_string( num_entities() ) );
		}

		if (( linestream >> PDBnum_string ).fail() ) {
			throw utility::excn::EXCN_Msg_Exception(
				"Failed to read PDB id for a residue on line " + utility::to_string( line_count ) +
				" of the EntityCorrespondence file:\n" + line );
		}
		/* resfile code */
		bool condition;
#ifndef _WIN32
		condition = std::isalpha(*PDBnum_string.rbegin());
#else
		condition = isalpha(*PDBnum_string.rbegin());
#endif
			if (condition) {
			/// The last character is the insertion code.
			for ( Size ch = 0; ch < PDBnum_string.length()-1; ++ch ) {
				if ( ! isdigit( PDBnum_string[ ch ] ) && ( ch != 0 || PDBnum_string[ ch ] != '-'  )) {
					throw utility::excn::EXCN_Msg_Exception(
						"Failed to read PDB residue id on line " + utility::to_string( line_count ) +
						" of the EntityCorrespondence file:\n" + line +"\nCharacter" +
						utility::to_string( ch + 1 ) + " of '" + PDBnum_string + "' is not a digit");
				}
			}
			PDBnum = atoi( PDBnum_string.substr( 0, PDBnum_string.length() - 1 ).c_str() );
			icode = *PDBnum_string.rbegin();
		} else { // no insertion code
			for ( Size ch = 0; ch < PDBnum_string.length(); ++ch ) {
				if ( ! isdigit( PDBnum_string[ ch ] ) && ( ch != 0 || PDBnum_string[ ch ] != '-'  )) {
					throw utility::excn::EXCN_Msg_Exception(
						"Failed to read PDB residue id on line " + utility::to_string( line_count ) +
						" of the EntityCorrespondence file:\n" + line +"\nCharacter" +
						utility::to_string( ch + 1 ) + " of '" + PDBnum_string + "' is not a digit");
				}
			}
			PDBnum = atoi( PDBnum_string.c_str() );
		}
		if ((linestream >> chain ).fail() ) {
			throw utility::excn::EXCN_Msg_Exception(
				"Failed to read a chain id for a PDB residue on line " + utility::to_string( line_count ) +
				" of the EntityCorrespondence file:\n" + line );
		}
		if (chain == '_') chain = ' ';
		residue_index = pdb_pose_map_->find( chain, PDBnum, icode );
		/* end resfile code */

		if ( residue_index == 0 ) {
			throw utility::excn::EXCN_Msg_Exception(
				"Residue ID read on line " + utility::to_string( line_count ) +
				" of the EntityCorrespondence file is not present in the pose: " +
				PDBnum_string + " ch: " + chain + " vs pose.total_residue()= " + utility::to_string( num_residues() ) );
		}
		add_resid_to_entity_list( entity_id, residue_index );
	}

	/*if ( count_correspondences == 0 ) {
		throw utility::excn::EXCN_Msg_Exception(
			"Failed to find any correspondences while reading from EntityCorrespondence file" );
			}*/
}

void
EntityCorrespondence::add_resid_to_entity_list(
	Size EntityID,
	Size ResID
)
{
	bounds_check_residue( funcnames_[ func_add_resid_to_entity_list ], ResID );
	bounds_check_entity( funcnames_[ func_add_resid_to_entity_list ], EntityID );
	if ( resid_2_entity_id_[ ResID  ] != 0 ) {
		if ( resid_2_entity_id_[ ResID  ] == EntityID ) return;
		entity_id_2_resids_[ resid_2_entity_id_[ ResID  ] ].remove( ResID );
	}
	entity_id_2_resids_[ EntityID ].push_back( ResID );
	resid_2_entity_id_[ ResID ] = EntityID;
}


EntityCorrespondence::Size
EntityCorrespondence::num_entities() const
{
	return entity_id_2_resids_.size();
}


EntityCorrespondence::Size
EntityCorrespondence::num_residues() const
{
	return pose_->total_residue();
}

EntityCorrespondence::Size
EntityCorrespondence::entity_for_residue( Size resid ) const
{
	bounds_check_residue( funcnames_[ func_entity_for_residue ], resid );
	return resid_2_entity_id_[ resid ];
}

EntityCorrespondence::Size
EntityCorrespondence::n_residues_for_entity( Size entity_id ) const
{
	bounds_check_entity( funcnames_[ func_n_residues_for_entity ], entity_id );
	return entity_id_2_resids_[ entity_id ].size();
}

EntityCorrespondence::ResIDListConstIter
EntityCorrespondence::residues_for_entity_begin( Size entity_id ) const
{
	bounds_check_entity( funcnames_[ func_residues_for_entity_begin ], entity_id );
	return entity_id_2_resids_[ entity_id ].begin();

}

EntityCorrespondence::ResIDListConstIter
EntityCorrespondence::residues_for_entity_end( Size entity_id ) const
{
	bounds_check_entity( funcnames_[ func_residues_for_entity_end ], entity_id );
	return entity_id_2_resids_[ entity_id ].end();

}

void EntityCorrespondence::bounds_check_entity(
	std::string const & funcname,
	Size entity_id
) const
{
	if ( entity_id > num_entities() ) {
		throw utility::excn::EXCN_Msg_Exception( "EntityCorrespondence::" + funcname + " "
			"attempted to access information for entity " + utility::to_string( entity_id ) +
			" (only " + utility::to_string( num_entities() ) + " entities exist)" );
	}
}

void EntityCorrespondence::bounds_check_residue(
	std::string const & funcname,
	Size resid
) const
{
	if ( ! pose_ ) {
		throw utility::excn::EXCN_Msg_Exception( "EntityCorrespondence::" + funcname + " "
			"invoked before pose_ pointer was set." );
	}
	if ( resid > num_residues() ) {
		throw utility::excn::EXCN_Msg_Exception( "EntityCorrespondence::" + funcname + " "
			"attempted to access information for residue " + utility::to_string( resid ) +
			" (only " + utility::to_string( num_residues() ) + " residues exist)" );
	}
}


}
}
