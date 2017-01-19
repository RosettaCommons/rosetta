// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/import_pose.hh
///
/// @brief  various functions to construct Pose object(s) from PDB(s)
/// @author Sergey Lyskov

#ifndef INCLUDED_core_import_pose_import_pose_hh
#define INCLUDED_core_import_pose_import_pose_hh

// Package headers
#include <core/import_pose/import_pose_options.fwd.hh>

// C++ headers
#include <iosfwd>

// Utility headers
#include <basic/Tracer.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

class CifFile;
class CifParser;
typedef utility::pointer::shared_ptr< CifFile > CifFileOP;
typedef utility::pointer::shared_ptr< CifParser > CifParserOP;


namespace core {
namespace import_pose {

enum FileType{
	PDB_file,
	CIF_file,
	Unknown_file
};

std::ostream & operator<<( std::ostream & stream, FileType type );

typedef std::string String;

void
read_all_poses(
	utility::vector1< std::string > const & filenames,
	utility::vector1< core::pose::PoseOP > &  poses
);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	io::StructFileRepCOP fd,
	bool read_fold_tree = false
);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);


/// @brief Returns a PoseOP object from the Pose created from input
/// PDB  <filename>
/// @note: in PyRosetta, this will return a Pose object
///
/// example(s):
///     pose = pose_from_file("YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
///     make_pose_from_sequence
///     pose_from_rcsb
///     pose_from_sequence
pose::PoseOP pose_from_file(
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Returns a PoseOP object from the Pose created from input
/// PDB  <filename>, taking a set of custom ImportPoseOptions parameters.
/// @note: in PyRosetta, this will return a Pose object
///
/// example(s):
///     pose = pose_from_file("YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
///     make_pose_from_sequence
///     pose_from_rcsb
///     pose_from_sequence
pose::PoseOP
pose_from_file(
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
);
/// @brief Returns a PoseOP object from the Pose created by reading the input
/// PDB  <filename>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
pose::PoseOP pose_from_file(
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Determine what file type is passed to function
/// there should only be one function that calls this, pose_from_file
/// and only calls it when the filetype is unknown
FileType
determine_file_type( std::string const &contents_of_file);

void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType file_type = Unknown_file
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, uses the FA_STANDARD ResidueTypeSet (fullatom) by default
/// @note: will use centroid if in::file::centroid_input is true
///
/// example(s):
///     pose_from_file(pose,"YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
void
pose_from_file(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_file(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>  using
/// ResidueTypeSet  <residue_set>
void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type = Unknown_file
);

utility::vector1< core::pose::Pose >
poses_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::Pose >
poses_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree,
	FileType type = Unknown_file
);


// FA_STANDARD residue set

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>
/// using the FA_STANDARD ResidueTypeSet (fullatom)
void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	std::string const & filename = ""
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	ImportPoseOptions const & options,
	std::string const & filename = ""
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options,
	std::string const & filename
);

void pose_from_pdb_stream(
	pose::Pose & pose,
	std::istream & pdb_stream,
	std::string const & filename,
	ImportPoseOptions const & options
);

// uses the CENTROID residue_set

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>  using the CENTROID ResidueTypeSet (centroid)
void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false
);

/// @brief Create pose object, using given StructFileRep object.
/// If PDB cleanin specified - it will be applied first.
/// Constructs a ImportPoseOptions object from the command line
void build_pose(
	io::StructFileRepOP fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
);

/// @brief Create pose object, using given StructFileRep object.
/// If PDB cleanin specified - it will be applied first
void build_pose(
	io::StructFileRepOP fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

/// @brief Create pose object, using given StructFileRep object.
/// No PDB cleanin will be appliend.
void build_pose_as_is(
	io::StructFileRepOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

void build_pose_as_is2(
	io::StructFileRepCOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	ImportPoseOptions const & options
);

} // namespace import_pose
} // namespace core

#endif // INCLUDED_core_import_pose_import_pose_HH
