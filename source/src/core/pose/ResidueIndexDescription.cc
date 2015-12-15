// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/pose/ResidueIndexDescription.cc
/// @brief  Two classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pose/ResidueIndexDescription.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>

namespace core {
namespace pose {

ResidueIndexDescription::ResidueIndexDescription() :
	unassigned_( true ),
	pose_numbered_( false ),
	pose_index_( 0 ),
	chain_( ' ' ),
	resindex_( 0 ),
	insertion_code_( ' ' )
{}

ResidueIndexDescription::ResidueIndexDescription(
	core::Size pose_index
) :
	unassigned_( false ),
	pose_numbered_( true ),
	pose_index_( pose_index ),
	chain_( ' ' ),
	resindex_( 0 ),
	insertion_code_( ' ' )
{}

ResidueIndexDescription::ResidueIndexDescription(
	char chain,
	int  resindex,
	char insertion_code
) :
	unassigned_( false ),
	pose_numbered_( false ),
	pose_index_ ( 0 ),
	chain_( chain ),
	resindex_( resindex ),
	insertion_code_( insertion_code )
{}

ResidueIndexDescription::~ResidueIndexDescription() {}

core::Size
ResidueIndexDescription::resolve_index(
	core::pose::Pose const & pose
) const
{
	// return bogus index, but do not trip an error.
	if ( unassigned_ ) { return 0; }
	if ( pose_numbered_ ) {
		if ( pose_index_ > pose.total_residue() ) {
			throw utility::excn::EXCN_Msg_Exception( "Residue index description exceeds the number of residues in the Pose: pose_index_ = " +
				utility::to_string( pose_index_ ) + " vs pose.total_residue() " + utility::to_string( pose.total_residue() ));
		}
		return pose_index_;
	} else {
		// THIS HAS NOT BEEN TESTED OR CAREFULLY CONSIDERED
		// ONLY A SKETCH OF CODE THAT COULD BE PUT HERE
		core::Size resid = pose.pdb_info()->pdb2pose().find( chain_, resindex_, insertion_code_ );
		if ( resid == 0 ) {
			// Dealing with the inability to return a null char constant by having two error messages.  LAAAAAAME.
			if ( insertion_code() != ' ' ) {
				throw utility::excn::EXCN_Msg_Exception( "Unable to find PDB residue " + utility::to_string( resindex_ ) + insertion_code_ + ' ' + "on chain " + chain_ + " in input Pose" );
			}
			throw utility::excn::EXCN_Msg_Exception( "Unable to find PDB residue " + utility::to_string( resindex_ ) + insertion_code_ + "on chain " + chain_ + " in input Pose" );
		}
		return resid;
	}
}

ResidueIndexDescriptionFromFile::ResidueIndexDescriptionFromFile() :
	ResidueIndexDescription(),
	linenum_( 0 )
{}

ResidueIndexDescriptionFromFile::ResidueIndexDescriptionFromFile(
	std::string fname,
	core::Size linenum,
	core::Size pose_index
) :
	ResidueIndexDescription( pose_index ),
	fname_( fname ),
	linenum_( linenum )
{}

ResidueIndexDescriptionFromFile::ResidueIndexDescriptionFromFile(
	std::string fname,
	core::Size linenum,
	char chain,
	int  resindex,
	char insertion_code
) :
	ResidueIndexDescription( chain, resindex, insertion_code ),
	fname_( fname ),
	linenum_( linenum )
{}

core::Size
ResidueIndexDescriptionFromFile::resolve_index(
	core::pose::Pose const & pose
) const
{
	// return bogus index, but do not trip an error.
	if ( unassigned() ) { return 0; }

	if ( pose_numbered() ) {
		if ( pose_index() > pose.total_residue() ) {
			throw utility::excn::EXCN_Msg_Exception( "Residue index description given on line " + utility::to_string( linenum_ ) +
				" of the file named " + fname_ + " exceeds the number of residues in the Pose: pose_index_ = " +
				utility::to_string( pose_index() ) + " vs pose.total_residue() " + utility::to_string( pose.total_residue() ));
		}
		return pose_index();
	} else {
		// THIS HAS NOT BEEN TESTED OR CAREFULLY CONSIDERED
		// ONLY A SKETCH OF CODE THAT COULD BE PUT HERE
		core::Size resid = pose.pdb_info()->pdb2pose().find( chain(), resindex(), insertion_code() );
		if ( resid == 0 ) {
			// Dealing with the inability to return a null char constant by having two error messages.  LAAAAAAME.
			if ( insertion_code() != ' ' ) {
				throw utility::excn::EXCN_Msg_Exception( "Unable to find PDB residue " + utility::to_string( resindex() ) + insertion_code() + ' ' + "on chain " + chain() + " in input Pose" );
			}
			throw utility::excn::EXCN_Msg_Exception( "Unable to find PDB residue " + utility::to_string( resindex() ) + insertion_code() + "on chain " + chain() + " in input Pose" );
		}
		return resid;
	}
}

// stupid visual studio has a different definition of std::isalpha than the C++ standard
// so I have to write my own version of isalpha...
bool
character_is_USA_letter( char c )
{
	return ( 'A' <= c && c <= 'Z' ) || ( 'a' <= c && c <= 'z' );
}

void
parse_PDBnum_icode(
	std::string const & token,
	std::string const & fname,
	Size const lineno,
	int & PDBnum,
	char & icode
)
{
	std::istringstream PDBnum_s;
	if ( character_is_USA_letter( *token.rbegin() ) ) {
		PDBnum_s.str(token.substr(0, token.length() - 1));
		icode = *token.rbegin();
	} else {
		PDBnum_s.str(token);
		icode = ' ';
	}

	char remaining;
	if ( !(PDBnum_s >> PDBnum) || PDBnum_s.get(remaining) ) {
		std::stringstream err_msg;
		err_msg
			<< "Problem in parse_PDBNum_icode. On line " << lineno
			<< " of the file '" << fname << "', "
			<< "the token '" << token << "' "
			<< "is not a valid <PDBNUM>[<ICODE>] identifier.";
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}

}

}
}
