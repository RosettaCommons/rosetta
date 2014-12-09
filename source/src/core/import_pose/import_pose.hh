// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <iostream>
// AUTO-REMOVED #include <string>

// Utility headers
#include <basic/Tracer.fwd.hh>
// AUTO-REMOVED #include <utility/io/ozstream.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.fwd.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {

/// @brief special Tracer instance acting as special param for all traced_dump_pdb functions
/// extern basic::Tracer TR_dump_pdb_dummy;

void
read_all_poses(
  const utility::vector1<std::string>& filenames,
	utility::vector1<core::pose::Pose>* poses);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	io::pdb::FileData const & fd,
	bool read_fold_tree = false
);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);

/* Undefined, commenting out to fix PyRosetta build  void
read_pdbinfo_labels(
  std::string const & s,
  pose::Pose & pose
); */

/// @brief Returns a PoseOP object from the Pose created from input
/// PDB  <filename>
/// @note: in PyRosetta, this will return a Pose object
///
/// example(s):
///     pose = pose_from_pdb("YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
///     make_pose_from_sequence
///     pose_from_rcsb
///     pose_from_sequence
pose::PoseOP pose_from_pdb(
	std::string const & filename,
	bool read_fold_tree = false
);

/// @brief Returns a PoseOP object from the Pose created by reading the input
/// PDB  <filename>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
pose::PoseOP pose_from_pdb(
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false
);

void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, uses the FA_STANDARD ResidueTypeSet (fullatom) by default
/// @note: will use centroid if in::file::centroid_input is true
///
/// example(s):
///     pose_from_pdb(pose,"YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
void
pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false
);

void
pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>  using
/// ResidueTypeSet  <residue_set>
void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false
);

void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
  chemical::ResidueTypeSet const & residue_set,
  utility::vector1< std::string > const & filenames,
  ImportPoseOptions const & options,
  bool read_fold_tree
);

utility::vector1< core::pose::Pose >
poses_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false
);

utility::vector1< core::pose::Pose >
poses_from_pdbs(
  chemical::ResidueTypeSet const & residue_set,
  utility::vector1< std::string > const & filenames,
  bool read_fold_tree
);

/* utility::vector1< core::pose::Pose >
poses_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
); */

/* utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false
); */

// FA_STANDARD residue set

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>
/// using the FA_STANDARD ResidueTypeSet (fullatom)
void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
	bool read_fold_tree = false
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

/// uses the CENTROID residue_set

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>  using the CENTROID ResidueTypeSet (centroid)
void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false
);

typedef std::string String;

void set_reasonable_fold_tree( core::pose::Pose & pose );

/// @brief Look for peptide bonds connected by >3A and replace them with a jump
// Undefinded commenting out to fix PyRosetta build void convert_missing_dens_to_jump( pose::Pose & pose );

/// @brief Create pose object, using given FileData object.
/// If PDB cleanin specified - it will be applied first.
/// Constructs a ImportPoseOptions object from the command line
void build_pose(
	io::pdb::FileData & fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
);

/// @brief Create pose object, using given FileData object.
/// If PDB cleanin specified - it will be applied first
void build_pose(
	io::pdb::FileData & fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

/// @brief Create pose object, using given FileData object.
/// No PDB cleanin will be appliend.
void build_pose_as_is(
	io::pdb::FileData & fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

} // namespace import_pose
} // namespace core

#endif // INCLUDED_core_import_pose_import_pose_HH
