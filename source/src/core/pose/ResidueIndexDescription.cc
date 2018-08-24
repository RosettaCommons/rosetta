// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pose/ResidueIndexDescription.cc
/// @brief  Classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/pose/ResidueIndexDescription.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>

namespace core {
namespace pose {

std::string RID_Source::source_string() const {
	return "";
}

std::string
RID_FileSource::source_string() const {
	return  "line " + utility::to_string( linenum_ ) + " of the file named " + fname_;
}

ResidueIndexDescription::~ResidueIndexDescription() = default;

std::string
ResidueIndexDescription::source_string() const {
	if ( source_ == nullptr ) {
		return "";
	} else {
		return "from " + source_->source_string() + " ";
	}
}

ResidueIndexDescriptionPoseNum::~ResidueIndexDescriptionPoseNum() = default;

core::Size
ResidueIndexDescriptionPoseNum::resolve_index(
	core::pose::Pose const & pose
) const
{
	if ( pose_index_ > pose.size() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Residue index description " + source_string() + "exceeds the number of residues in the Pose: pose_index_ = " +
			utility::to_string( pose_index_ ) + " vs pose.size() " + utility::to_string( pose.size() ));
	}
	return pose_index_;
}

ResidueIndexDescriptionPDB::~ResidueIndexDescriptionPDB() = default;

core::Size
ResidueIndexDescriptionPDB::resolve_index(
	core::pose::Pose const & pose
) const
{
	// THIS HAS NOT BEEN TESTED OR CAREFULLY CONSIDERED
	// ONLY A SKETCH OF CODE THAT COULD BE PUT HERE
	core::Size resid = pose.pdb_info()->pdb2pose().find( chain_, resindex_, insertion_code_ );
	if ( resid == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "PDB residue " + source_string() + utility::to_string( resindex_ ) + insertion_code_ + ' ' + "on chain " + chain_ + " " + source_string() + "is not present in the input Pose." );
	}
	return resid;
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
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg.str() );
	}

}

}
}
