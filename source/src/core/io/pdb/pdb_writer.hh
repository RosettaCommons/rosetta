// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/pdb/pdb_writer.hh
/// @brief Function(s) for pdb writing
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team


#ifndef INCLUDED_core_io_pdb_pdb_writer_hh
#define INCLUDED_core_io_pdb_pdb_writer_hh

// Unit headers
#include <core/io/pdb/Record.hh>

// Package headers
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>

// Project header
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility header
#include <utility/io/ozstream.fwd.hh>


namespace core {
namespace io {
namespace pdb {


/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
///  Use default StructFileRepOptions
bool
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name
);

/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
bool
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	core::io::StructFileRepOptionsCOP options
);


/// @brief Writes  <pose>  to a given stream in PDB file format
///  Use default StructFileRepOptions
/// If a non-empty string is given for the out_fname variable, then
/// the Pose energy table at the bottom of the PDB will be labeled
/// with this string.
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & out_fname = ""
);

/// @brief Writes  <pose>  to a given stream in PDB file format
/// If a non-empty string is given for the out_fname variable, then
/// the Pose energy table at the bottom of the PDB will be labeled
/// with this string.
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options,
	std::string const & out_fname = ""
);




/// @brief This version takes an AtomID mask.
/// @details Used by Will's motif hash stuff, I think.
void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	id::AtomID_Mask const & mask,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);



void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	utility::vector1< core::Size > const & residue_indices,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


/// @brief Writes poses to a single PDB file, returns false if an error occurs
///  Use default StructFileRepOptions
bool
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::string const & file_name
);

/// @brief Writes poses to a single PDB file, returns false if an error occurs
bool
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::string const & file_name,
	core::io::StructFileRepOptionsCOP options
);

/// @brief Writes poses to a given stream in PDB file format
///  Use default StructFileRepOptions
void
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::ostream & out
);

/// @brief Writes poses to a given stream in PDB file format
void
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options
);


/// @brief Adds a pose to a multimodel pdb file (or creates the file)
///  Use default StructFileRepOptions
bool
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	std::string const & model_tag,
	bool clear_existing_structures = false
);

/// @brief Adds a pose to a multimodel pdb file (or creates the file)
bool
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	std::string const & model_tag,
	core::io::StructFileRepOptionsCOP options,
	bool clear_existing_structures = false
);

/// @brief Adds a pose to a multimodel pdb file
///  Use default StructFileRepOptions
void
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & model_tag
);

/// @brief Adds a pose to a multimodel pdb file
void
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & model_tag,
	core::io::StructFileRepOptionsCOP options
);




/// @brief PyRosetta-Compatible; Writes pdb data for the given residue, beginning from the given atom number
std::string
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions ),
	core::Size start_atom_number = 1
);

/// @brief PyRosetta-Compatible; Writes pdb data for the given residue, incrementing atom_number counter
std::string
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::Size & atom_number,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);


/// @brief Writes pdb data for the given residue, incrementing atom_number counter
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::Size & atom_number,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options=core::io::StructFileRepOptionsCOP( new core::io::StructFileRepOptions )
);

/// @brief Writes pdb data for the given residue, beginning from the given atom number
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	std::ostream & out,
	core::Size start_atom_number = 1,
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





///////////////////////////////////////////////////////////////
//////// LEGACY JD2-layer compatability (DO NOT USE) //////////
///////////////////////////////////////////////////////////////


/// @brief Writes a pose to a given stream in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] jd2_job_data Additional data to append to the PDB file data (from the job).
/// @param[out] out The output stream that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
///
/// NOTE: The two booleans are now specified in the StructFileRepOptions. Job data is no longer used in JD3. Do not use this function for general pdb writing.  Use dump_file.
///
/// @author Vikram K. Mulligan (vmullig@uw.edu).
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & jd2_job_data,
	utility::io::ozstream & out,
	std::string const &filename=""
);

/// @brief Writes a pose to a given string in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] jd2_job_data Additional data to append to the PDB file data from the Job.
/// @param[out] out The output string that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
///
/// NOTE: The two booleans are now specified in the StructFileRepOptions. Job data is no longer used in JD3. Do not use this function for general pdb writing.  Use dump_file.
///
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & jd2_job_data,
	std::string & out,
	std::string const &filename=""
);
} //pdb
} //io
} //core


#endif //core/io/pdb/pdb_writer_hh

