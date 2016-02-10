// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/io/pdb/pdb_writer.hh
/// @brief Function(s) for pdb writing
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team


#ifndef INCLUDED_core_io_pdb_pdb_writer_hh
#define INCLUDED_core_io_pdb_pdb_writer_hh

#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/Field.fwd.hh>

// Package header
#include <core/io/StructFileRep.hh>

// Utility header
#include <utility/io/ozstream.fwd.hh>


namespace core {
namespace io {
namespace pdb {

/// @brief Writes a pose to a given stream in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] extra_data Additional data to append to the PDB file data.
/// @param[in] add_score_data Grab additional score data from the pose?
/// @param[in] add_extra_score_data Grab still more score data from the pose?
/// @param[out] out The output stream that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &extra_data,
	bool const add_score_data,
	bool const add_extra_score_data,
	utility::io::ozstream & out,
	std::string const &filename="",
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief This version takes an AtomID mask.
/// @details Used by Will's motif hash stuff, I think.
void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	id::AtomID_Mask const & mask,
	std::string const &tag,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief Writes a pose to a given string in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] extra_data Additional data to append to the PDB file data.
/// @param[in] add_score_data Grab additional score data from the pose?
/// @param[in] add_extra_score_data Grab still more score data from the pose?
/// @param[out] out The output string that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &extra_data,
	bool const add_score_data,
	bool const add_extra_score_data,
	std::string & out,
	std::string const &filename="",
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief Writes  <pose>  to a given stream in PDB file format
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & tag="",
	bool write_fold_tree = false,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
bool
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	std::string const & tag="",
	bool write_fold_tree = false,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	utility::vector1< core::Size > const & residue_indices,
	std::string const & tag="",
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


/// @brief Writes pdb data for the given residue, incrementing atom_number counter
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	Size & atom_number,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief Writes pdb data for the given residue, beginning from the given atom number
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	std::ostream & out,
	Size start_atom_number = 1,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


/// @brief Create a full .pdb as a string given a StructFileRep object.
std::string create_pdb_contents_from_sfr(
	StructFileRep const & sfr,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief Create a formatted line for .pdb output.
std::string create_pdb_line_from_record( Record const & record );

/// @brief Create vector of records from given StructFileRep object.
std::vector< Record > create_records_from_sfr(
	StructFileRep const & sfr,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

} //pdb
} //io
} //core


#endif //core/io/pdb/pdb_writer_hh

