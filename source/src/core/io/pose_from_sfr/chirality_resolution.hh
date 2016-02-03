// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pose_from_sfr/chirality_resolution.hh
/// @brief  Various utilities to accomodated PDB input issues.
/// @author Sergey Lyskov, Rhiju Das, Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_io_pose_from_sfr_chirality_resolution_hh
#define INCLUDED_core_io_pose_from_sfr_chirality_resolution_hh

// Unit headers
#include <core/io/pose_from_sfr/chirality_resolution.fwd.hh> // For typedefs

// Package headers
#include <core/io/StructFileRep.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>

//External headers

#include <boost/bimap.hpp>

// C++ headers
#include <iostream>
#include <map>
#include <string>


// Some  useful functions below for fixing up H5'<-->H5'' ambiguities based on geometric
//  comparisons to ideal coordinates.


namespace core {
namespace io {
namespace pose_from_sfr {

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief due to differences in different crystallography/NMR/modeling packages, labeling of sister atoms
///  (like OP1 <--> OP2, or H41 <--> H42) in PDBs is totally wacky. This is an attempt to regularize...
///  and it can actually make a difference since sometimes partial charges on sister hydrogens
///  can be different. Right now only set up for nucleic acids, but could probably generalize.
void check_and_correct_sister_atoms( core::conformation::ResidueOP & rsd );

void check_and_correct_sister_atom_based_on_chirality( core::conformation::ResidueOP & rsd,
	std::string const & sister1_name,
	std::string const & sister2_name,
	std::string const & parent_name,
	std::string const & outer_ref_name );

void
check_and_correct_sister_atom_based_on_outgroup( core::conformation::ResidueOP & rsd,
	std::string const & sister1_name,
	std::string const & sister2_name,
	std::string const & outgroup_name );


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


/// @brief Get theshold distance below which two atoms are considered bonded. (1.2*covalent)
core::Real
bonding_distance_threshold( std::string element1, std::string element2 );

/// @brief Scoring scheme for the heuristic PDB renaming
core::Real
score_mapping( NameBimap const & mapping,
	ResidueInformation const & rinfo,
	chemical::ResidueType const & rsd_type );

/// @brief Attempt to use element identity and connectivity to map atom names from the rinfo object onto the rsd_type object names.
void
remap_names_on_geometry( NameBimap & mapping,
	ResidueInformation const & rinfo,
	chemical::ResidueType const & rsd_type);


} // namespace pose_from_sfr
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pose_from_sfr_chirality_fixup_HH
