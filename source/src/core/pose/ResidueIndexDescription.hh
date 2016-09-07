// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pose/ResidueIndexDescription.hh
/// @brief  Two classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pose_ResidueIndexDescription_HH
#define INCLUDED_core_pose_ResidueIndexDescription_HH

// Unit headers
#include <core/pose/ResidueIndexDescription.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace core {
namespace pose {

/// @brief a class which can represent one of many ways in which to describe a
/// particular residue in a pose, and can, when given a pose, find its index.
/// The object should be constructed with all its needed parameters, but, one
/// instance may be copied from another.
class ResidueIndexDescription : public utility::pointer::ReferenceCount
{
public:
	ResidueIndexDescription();
	ResidueIndexDescription( Size pose_index );
	ResidueIndexDescription(
		char chain,
		int resindex,
		char insertion_code = ' ' // space character represents no insertion code
	);
	~ResidueIndexDescription() override;

	/// @breif Turn the internal data into a pose index (in the range from 1 to total residue)
	/// or throw an exception if the index cannot be resolved.
	virtual Size resolve_index( core::pose::Pose const & p ) const;

	bool unassigned() const { return unassigned_; }
	bool pose_numbered() const { return pose_numbered_; }
	Size  pose_index() const { return pose_index_; }
	char chain() const { return chain_; }
	int resindex() const { return resindex_; }
	char insertion_code() const { return insertion_code_; }

private:
	bool unassigned_;
	bool pose_numbered_;
	Size  pose_index_;
	char chain_;   // chain character
	int resindex_;
	char insertion_code_;
};

/// @brief %ResidueIndexDescriptionFromFile differs from its parent only
/// in its error reporting when a residue index resolution fails: it
/// prints the file name and the line number that the residue index data
/// came from.
class ResidueIndexDescriptionFromFile : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionFromFile();

	ResidueIndexDescriptionFromFile(
		std::string fname,
		Size linenum,
		Size pose_index
	);

	ResidueIndexDescriptionFromFile(
		std::string fname,
		Size linenum,
		char chain,
		int resindex,
		char insertion_code = ' ' // space character represents no insertion code
	);

	/// @breif Turn the internal data into a pose index (in the range from 1 to total residue)
	/// or throw an exception if the index cannot be resolved; the exception thrown by this
	/// class if the resolution fails includes a statement of the file name and the line
	/// number in that file.
	Size resolve_index( core::pose::Pose const & p ) const override;

	std::string const & fname() const { return fname_; }
	Size linenum() const { return linenum_; }

private:
	std::string fname_; // the file the residue index descriptor was declared in; for error reporting
	Size linenum_; // the line number the residue index descriptor was declared in; for error reporting
};

/// @brief Take the string "token" and try to interpret it as a PDB identifier in the form of an
/// integer as well as an optional insertion code.  For example the string "25A" would be
/// interpretted as the residue 25 with the insertion code "A."  Throws an exception if
/// the input string is misformatted.
void
parse_PDBnum_icode(
	std::string const & token,
	std::string const & fname,
	Size const lineno,
	int & PDBnum,
	char & icode
);


}
}

#endif
