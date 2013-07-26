// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.hh
///
/// @brief
/// @author Sergey Lyskov

#ifndef INCLUDED_core_io_pdb_file_data_fixup_hh
#define INCLUDED_core_io_pdb_file_data_fixup_hh


// Unit headers
#include <core/io/pdb/file_data_fixup.hh>

// Package headers
#include <core/io/pdb/file_data_options.fwd.hh>
#include <core/io/pdb/file_data.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>

// C++ headers
#include <iostream>
#include <map>
#include <string>


////////////////////////////////////////////////////////////////////////////////////////////
// useful utility scripts for fixing up residue and atom names.
//
// Originally implemented for nucleic acids, to permit backward-compatibility with
// old Rosetta atom & residue names.
//
// Some  useful functions below for fixing up H5'<-->H5'' ambiguities based on geometric
//  comparisons to ideal coordinates.
//
//

namespace core {
namespace io {
namespace pdb {

std::string
convert_res_name( std::string const & name );

std::string
convert_atom_name( std::string const & res_name, std::string atom_name );

void convert_nucleic_acid_residue_info_to_standard( 	utility::vector1< ResidueInformation > & rinfo );

bool is_potential_old_DNA( std::string const & res_name );

bool is_old_RNA( std::string const & res_name );

bool is_NA( std::string const & res_name );

bool missing_O2prime( utility::vector1< AtomInformation > const & atoms );

// @brief  This is a pretty good framework and could allow for other crazy nucleic acid atom name schemes.
void convert_nucleic_acid_atom_names_to_standard( ResidueInformation & rinfo  );

void convert_nucleic_acid_atom_name_to_standard( AtomInformation & atom_info );

//////////////////////////////////////////////////////////////////////////////////////////
// @brief due to differences in different crystallography/NMR/modeling packages, labeling of sister atoms
//  (like OP1 <--> OP2, or H41 <--> H42) in PDBs is totally wacky. This is an attempt to regularize...
//  and it can actually make a difference since sometimes partial charges on sister hydrogens
//  can be different. Right now only set up for nucleic acids, but could probably generalize.
void check_and_correct_sister_atoms( core::conformation::ResidueOP & rsd );

void check_and_correct_sister_atom_based_on_chirality( core::conformation::ResidueOP & rsd,
																											 std::string const sister1_name,
																											 std::string const sister2_name,
																											 std::string const parent_name,
																											 std::string const outer_ref_name );

void
check_and_correct_sister_atom_based_on_outgroup( core::conformation::ResidueOP & rsd,
																								 std::string const sister1_name,
																								 std::string const sister2_name,
																								 std::string const outgroup_name );


void flip_atom_xyz( core::conformation::ResidueOP & rsd,
										std::string const & sister1_name,
										std::string const & sister2_name );

int sgn( Real const & x );

int
get_chirality_sign(  Vector const & xyz_sister1,
										 Vector const & xyz_sister2,
										 Vector const & xyz_parent,
										 Vector const & xyz_outer_ref );

int
get_closest_sister(  Vector const & xyz_sister1,
										 Vector const & xyz_sister2,
										 Vector const & xyz_outgroup );


} // namespace pdb
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pdb_file_data_HH
