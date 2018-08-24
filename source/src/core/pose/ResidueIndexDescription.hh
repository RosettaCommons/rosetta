// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pose/ResidueIndexDescription.hh
/// @brief  Classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

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

/// @brief A class to hold information about where a ResidueIndexDescriptor comes from,
/// to allow for better error messages
class RID_Source : public utility::pointer::ReferenceCount
{
public:
	virtual std::string source_string() const;
};


/// @brief Representation of a file input support A class to hold information about where a ResidueIndexDescriptor comes from,
/// to allow for better error messages
class RID_FileSource : public RID_Source
{
public:
	RID_FileSource(
		std::string const & fname,
		Size linenum
	):
		fname_(fname),
		linenum_(linenum)
	{}

	std::string const & fname() const { return fname_; }
	core::Size linenum() const { return linenum_; }

	std::string source_string() const override;

private:
	std::string fname_;
	Size linenum_ = 0;
};

/// @brief a class which can represent one of many ways in which to describe a
/// particular residue in a pose, and can, when given a pose, find its index.
/// The object should be constructed with all its needed parameters, but, one
/// instance may be copied from another.
class ResidueIndexDescription : public utility::pointer::ReferenceCount
{
public:
	ResidueIndexDescription(): source_(nullptr) {}
	ResidueIndexDescription(RID_SourceCOP source): source_(source) {}

	~ResidueIndexDescription() override;

	/// @brief Turn the internal data into a pose index (in the range from 1 to total residue)
	/// or throw an exception if the index cannot be resolved.
	virtual Size resolve_index( core::pose::Pose const & p ) const = 0;

	RID_SourceCOP get_source() const { return source_; }

protected:

	/// @brief Provide a formatted string representing the source of the information (if any)
	std::string source_string() const;

private:
	RID_SourceCOP source_;
};

/// @brief a class which represents a residue index as a literal, Rosetta/Pose numbered integer
class ResidueIndexDescriptionPoseNum : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionPoseNum( Size pose_index ):
		pose_index_(pose_index)
	{}
	ResidueIndexDescriptionPoseNum(RID_SourceCOP source, Size pose_index ):
		ResidueIndexDescription(source),
		pose_index_(pose_index)
	{}

	~ResidueIndexDescriptionPoseNum() override;

	Size pose_index() const { return pose_index_; }

	Size resolve_index( core::pose::Pose const & p ) const override;

private:
	Size  pose_index_;
};

/// @brief a class which represents a residue index as a PDB information (chain, resindex, insertion code)
class ResidueIndexDescriptionPDB : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionPDB(
		char chain,
		int resindex,
		char insertion_code = ' ' // space character represents no insertion code
	) :
		chain_(chain),
		resindex_(resindex),
		insertion_code_(insertion_code)
	{}
	ResidueIndexDescriptionPDB(
		RID_SourceCOP source,
		char chain,
		int resindex,
		char insertion_code = ' ' // space character represents no insertion code
	) :
		ResidueIndexDescription(source),
		chain_(chain),
		resindex_(resindex),
		insertion_code_(insertion_code)
	{}

	~ResidueIndexDescriptionPDB() override;

	char chain() const { return chain_; }
	int resindex() const { return resindex_; }
	char insertion_code() const { return insertion_code_; }

	Size resolve_index( core::pose::Pose const & p ) const override;

private:
	char chain_;   // chain character
	int resindex_;
	char insertion_code_;
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
