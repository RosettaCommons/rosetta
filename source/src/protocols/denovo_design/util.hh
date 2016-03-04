// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/util.hh
/// @brief util functions for building structures from components
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_util_hh
#define INCLUDED_protocols_denovo_design_util_hh

// Unit headers

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Package headers

// Core headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/Remarks.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <basic/datacache/DataMap.fwd.hh>

// C++ Headers
#include <set>
#include <list>
#include <map>

namespace protocols {
namespace denovo_design {

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & set1,
	std::set< core::Size > const & set2,
	bool const keep_chirality=false
);

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & res_set,
	bool const keep_chirality=false
);

/// @brief builds an extended-conformation poly-XXX pose
core::pose::PoseOP construct_dummy_pose( std::string const & restype_name );
core::pose::PoseOP construct_dummy_pose( std::string const & restype_name, core::Size const length );
core::pose::PoseOP construct_dummy_pose( core::chemical::ResidueType const & restype );
core::pose::PoseOP construct_dummy_pose( core::chemical::ResidueType const & restype, core::Size const length );

core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( basic::datacache::DataMap const & data, std::string const & name );

std::string
abego_str( utility::vector1< std::string > const & abego );

utility::vector1< std::string >
abego_vector( std::string const & ab );

// gets a remark line, pasting multiple lines together if necessary
std::string get_remark_line( core::io::Remarks::const_iterator & it_rem, core::io::Remarks::const_iterator const & end );

/// @brief adds a remark to a Remarks object, splitting it into multiple remarks if it is too long
void add_remark( core::io::Remarks & remarks, core::Size const num, std::string const & str_val );

// helper function to calculate stop of loop without overlap
core::Size
loop_stop_without_overlap( core::pose::Pose const & pose, core::Size stopres, core::Size const overlap );

// helper function to calculate stop residue of loop without overlap
core::Size
loop_start_without_overlap( core::pose::Pose const & pose, core::Size startres, core::Size const overlap );

core::kinematics::FoldTree
remove_jump_atoms( core::kinematics::FoldTree const & orig );

/// @brief given a residue, rebuilds all missing atoms
void rebuild_missing_atoms( core::pose::Pose & pose, core::Size const resi );

/// @brief helper function that looks for the given residue in a fold tree and returns the jump that controls its 6D-DoFs
int find_jump_rec(
	core::kinematics::FoldTree const & ft,
	int const residue );

/// @brief inserts the peptide edges to accomodate the new jump edge given
void insert_peptide_edges( core::kinematics::FoldTree & ft, core::kinematics::Edge const & jedge );

// @brief parses a string containing single integers and ranges. Returns a vector of all possible values.
utility::vector1< core::Size >
parse_length_string( std::string const & len_str );

/// @brief given a number 0 <= x < 1, calculate an integer M <= x <= N
/// NOTE THAT THIS FUNCTION MODIFIES THE PARAMETER
core::Size
extract_int( core::Real & num, core::Size const m, core::Size const n );

/// @brief copies rotamers from the pose "src" into the permutation "dest"
/// no backbone changes are made
/// if detect_disulf flag is on, disulfides will be re-detected
void
copy_rotamers( components::StructureData & dest, core::pose::Pose const & src );

/// @brief counts the beta-bulges in the given segment. This simply counts all abego A's in the segment
core::Size
count_bulges( components::StructureData const & perm, std::string const & segment );

/// @brief gets all strand pairings from a perm
/// @details if use_register_shift=0, the returned register shift is 99
std::string get_strandpairings(
	components::StructureData const & perm,
	bool const use_register_shift );

/// @brief dumps a pose into another pose as a new chain
void
add_chain_from_pose( core::pose::PoseCOP to_add, core::pose::PoseOP combined );

} // denovo_design
} // protocols

//////////////////////////////////////////////////////////////////////////
/// Output operators for std template classes                          ///
//////////////////////////////////////////////////////////////////////////
namespace std {

/// @brief outputs a set
std::ostream & operator<<( std::ostream & os, std::set< int > const & set );

/// @brief outputs a set
std::ostream & operator<<( std::ostream & os, std::set< core::Size > const & set );

/// @brief outputs a set
std::ostream & operator<<( std::ostream & os, std::set< std::string > const & set );

/// @brief outputs a list of sizes
std::ostream & operator<<( std::ostream & os, std::list< core::Size > const & list );

/// @brief outputs a list of strings
std::ostream & operator<<( std::ostream & os, std::list< std::string > const & list );

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< core::Size, core::Size > const & map );

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< char, core::Size > const & map );

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< std::string, core::Size > const & map );

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< std::pair< std::string, std::string >, core::Size > const & map );

/// @brief outputs a map
std::ostream & operator<<( std::ostream & os, std::map< std::string, core::Real > const & map );

/// @brief outputs a vector
std::ostream & operator<<( std::ostream & os, numeric::xyzVector< core::Real > const & vec );

} // std
#endif
