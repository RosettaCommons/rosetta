// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/build_pose_as_is.hh
/// @brief  Declarations for StructFileRep and related classes.
/// @author Sergey Lyskov

// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean the Resource Manager may not work properly.
// Instead, modify StructFileRepOptions to include the option.


#ifndef INCLUDED_core_io_pdb_file_data_hh
#define INCLUDED_core_io_pdb_file_data_hh


// Unit headers
#include <core/io/pdb/Field.fwd.hh>

// Package headers
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/io/izstream.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <map>
#include <string>


namespace core {
namespace io {
namespace pdb {

typedef std::string String;


/// @brief Builds a pose into  <pose>, without repacking or optimizing
/// hydrogens; using the full-atom ResidueTypeSet
void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename
);

/// @brief Builds a pose into  <pose>, without repacking or optimizing
/// hydrogens; using the full-atom ResidueTypeSet and a set of options.
void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename,
	StructFileReaderOptions const & pdr_options
);

void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
);

void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	StructFileReaderOptions const & pdr_options
);

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	std::istream & file_contents,
	StructFileReaderOptions const & pdr_options
);

} // namespace pdb
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_file_data_HH
