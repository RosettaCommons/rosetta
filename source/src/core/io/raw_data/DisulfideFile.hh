// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

/// @file   core/io/raw_data/DisulfideFile.hh
/// @brief  A simple wrapper for a Disulfide File suitable for the -fix_disulf option.
/// @author Spencer Bliven <blivens@u.washington.edu>

#ifndef INCLUDED_core_io_raw_data_DisulfideFile_hh
#define INCLUDED_core_io_raw_data_DisulfideFile_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <string>
// AUTO-REMOVED #include <utility>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {

/// @brief Parses and stores a disulfide file.
/// @details Initiallizing a DisulfideFile is a lightweight operation.
/// The heavy lifting occurs the first time disulfides() is called. This parses
/// the file and caches the resulting pairs of residues.
/// Subsequent calls to disulfides() are fast since they don't reparse the file
/// but merely reinterpret the results in terms of the specified.
///
/// @section format File Format
/// The disulfide file format is pretty flexible. It looks for three kinds of lines:
///  - '12	42' Lines with two integers are interpreted as a disulfide
///    bond, indexed by internal rosetta residue number, i.e. the 12th and 42nd
///    residues from the start of the pose.
///  - '12A	42A' If a single character is appended to the numbers, DisulfideFile
///    assumes that these refer to the pdb number and chain, i.e. the residues
///    of chain A numbered 12 and 42 in the pdb.
///  - 'SSBOND   1 CYS A   12    CYS A   42' The
///	   <a href="http://www.wwpdb.org/documentation/format32/sect6.html#SSBOND">
///    PDB format</a> for disulfide bond annotations.
///
/// All lines not matching these criteria are silently ignored. This implies
/// that whole pdb files can usually be used unaltered as disulfide files, since
/// the SSBOND entries are extracted and all else is ignored.
///
/// Question from rhiju, 2014 -- I don't see handling of '12A' or SSBOND cases in
///  DisulfideFile.cc -- is someone going to check that in?
///

class DisulfideFile {
private:
	/// @brief distinguish between PDB numbering and internal Rosetta numbering
	/// @details unknown_num should be avoided, as one of the other types can
	/// generally be infered from the existance of a chain specifier.
	enum NumberingSystem {
		pdb_num,
		rosetta_num,
		unknown_num
	};
	/// @brief represents a residue of either pdb or rosetta numbering.
	/// @details pdb numbers should include a chain character; rosetta numbers
	/// have chain==0
	typedef struct {
		core::Size n;
		char chain;
		NumberingSystem type;
	} ResNum;

public:
	DisulfideFile(std::string filename) :
		filename_(filename),
		up_to_date_(false),
		disulfides_() { }

	/// @brief Accessor for the filename
	inline std::string const& filename() const { return filename_; }

	/// @brief Get a list of disulfide bonds declared in the file
	void disulfides(
		utility::vector1< std::pair<core::Size,core::Size> > & disulfides ) const;

	/// @brief Get a list of disulfide bonds declared in the file
	///        (renumbered to rosetta numbering if necessary)
	void disulfides(
		utility::vector1< std::pair<core::Size,core::Size> > & disulfides,
		core::pose::Pose const& pose ) const;

	/// @brief Get a list of disulfide bonds declared in the file
	///        (renumbered to rosetta numbering if necessary)
    /// also manually set the disulfides in the conformation of the provided pose
    /// (this is a necessary workaround for dealing with multiple disulfide specification
    /// files in PyRosetta
    void read_in_and_set_disulfides(
        core::pose::Pose &pose );
private:
	/// @brief helper function to read in the file
	void parse_disulf_file() const;
	/// @brief convert a ResNum struct into a normal rosetta residue num
	core::Size resnum_to_rosetta_num(core::pose::Pose const& pose, ResNum const& resnum) const;
private:
	const std::string filename_;

	//Read from the file once, then cache data in disulfides_
	mutable bool up_to_date_; //flag that disulfides_ has read the file
	mutable utility::vector1< std::pair<ResNum,ResNum> > disulfides_;

}; //DisulfideFile

} //raw_data
} //io
} //core

#endif
