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
#include <iosfwd>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

/// @brief A class to hold information about where a ResidueIndexDescriptor comes from,
/// to allow for better error messages
class RID_Source : public utility::pointer::ReferenceCount
{
public:
	virtual std::string source_string() const;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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

#ifdef    SERIALIZATION
protected:
	RID_FileSource() = default;
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	/// If it can't be converted, throw an exception, unless no_error is true.
	/// If no_error is true, return 0 or a residue number outside the pose range.
	virtual Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const = 0;

	/// @brief Output a description of the residue.
	/// @details In general will give a form which can be added "in-line" (no line ending)
	virtual void show( std::ostream & ) const = 0;

	RID_SourceCOP get_source() const { return source_; }

protected:

	/// @brief Provide a formatted string representing the source of the information (if any)
	std::string source_string() const;

	/// @brief Convenience function for triggering error conditions
	core::Size do_error( bool no_error, std::string const & msg ) const;

private:
	RID_SourceCOP source_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

std::ostream & operator<<( std::ostream & out, ResidueIndexDescription const & rid );

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

	Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const override;

	void show( std::ostream & ) const override;

private:
	Size  pose_index_;

#ifdef    SERIALIZATION
protected:
	ResidueIndexDescriptionPoseNum() = default;
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief Convenience function for converting a Size (Pose numbered) into a ResidueIndexDescription
ResidueIndexDescriptionCOP
make_rid_posenum( core::Size resnum );

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

	Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const override;

	void show( std::ostream & ) const override;

private:
	char chain_;   // chain character
	int resindex_;
	char insertion_code_;

#ifdef    SERIALIZATION
protected:
	ResidueIndexDescriptionPDB() = default;
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief a class which represents a residue index as a reference-pose enabled information.
class ResidueIndexDescriptionRefPose : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionRefPose(
		std::string const & refpose_name,
		core::Size refpose_number,
		signed long refpose_offset
	):
		refpose_name_(refpose_name),
		refpose_number_(refpose_number),
		refpose_offset_(refpose_offset)
	{}

	ResidueIndexDescriptionRefPose(
		RID_SourceCOP source,
		std::string const & refpose_name,
		core::Size refpose_number,
		signed long refpose_offset
	):
		ResidueIndexDescription(source),
		refpose_name_(refpose_name),
		refpose_number_(refpose_number),
		refpose_offset_(refpose_offset)
	{}

	~ResidueIndexDescriptionRefPose() override;

	std::string const & refpose_name() const { return refpose_name_; }
	core::Size refpose_number() const { return refpose_number_; }
	signed long refpose_offset() const { return refpose_offset_; }

	Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const override;

	void show( std::ostream & ) const override;

private:
	std::string refpose_name_;
	core::Size refpose_number_;
	signed long refpose_offset_; //Must be signed!

#ifdef    SERIALIZATION
protected:
	ResidueIndexDescriptionRefPose() = default;
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief a class which represents the last residue in the Pose
class ResidueIndexDescriptionLastResidue : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionLastResidue()
	{}

	ResidueIndexDescriptionLastResidue(
		RID_SourceCOP source
	):
		ResidueIndexDescription(source)
	{}

	~ResidueIndexDescriptionLastResidue() override;

	Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const override;

	void show( std::ostream & ) const override;

private:
	// No data needed currently

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief a class which represents the last residue in the Pose
class ResidueIndexDescriptionChainEnd : public ResidueIndexDescription
{
public:
	ResidueIndexDescriptionChainEnd(char chain, bool chain_start = false):
		chain_start_( chain_start ),
		chain_letter_( chain )
	{}

	ResidueIndexDescriptionChainEnd(
		RID_SourceCOP source,
		char chain,
		bool chain_start = false
	):
		ResidueIndexDescription(source),
		chain_start_( chain_start ),
		chain_letter_( chain )
	{}

	ResidueIndexDescriptionChainEnd(core::Size chain, bool chain_start = false):
		chain_start_( chain_start ),
		chain_no_( chain )
	{}

	ResidueIndexDescriptionChainEnd(
		RID_SourceCOP source,
		core::Size chain,
		bool chain_start = false
	):
		ResidueIndexDescription(source),
		chain_start_( chain_start ),
		chain_no_( chain )
	{}

	~ResidueIndexDescriptionChainEnd() override;

	Size resolve_index( core::pose::Pose const & p, bool no_error = false ) const override;

	void show( std::ostream & ) const override;

private:
	/// @brief Should it be the front end of the chain we return?
	bool chain_start_ = false;
	/// @brief The chain number to use.
	/// If equal to -1, use chain letter instead.
	core::Size chain_no_ = -1;
	/// @brief If valid (non-null), identify chain by letter
	char chain_letter_ = 0x00; // This is NUL, not the character '0'

#ifdef    SERIALIZATION
protected:
	ResidueIndexDescriptionChainEnd() = default;
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_ResidueIndexDescription )
#endif // SERIALIZATION

#endif
