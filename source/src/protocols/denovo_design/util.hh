// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/string_util.hh>

// C++ Headers
#include <set>
#include <list>
#include <map>

namespace protocols {
namespace denovo_design {

typedef std::map< std::string, core::Size > StrandNameToNumberMap;

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
/// @details If keep_chirality is true, the D-amino acids are mutated to D-alanine.
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	core::select::residue_selector::ResidueSubset const & set1,
	core::select::residue_selector::ResidueSubset const & set2,
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
std::string
get_remark_line( core::io::Remarks::const_iterator & it_rem, core::io::Remarks::const_iterator const & end );

/// @brief adds a remark to a Remarks object, splitting it into multiple remarks if it is too long
void add_remark( core::io::Remarks & remarks, core::Size const num, std::string const & str_val );

// helper function to calculate stop of loop without overlap
core::Size
loop_stop_without_overlap( core::pose::Pose const & pose, core::Size stopres, core::Size const overlap );

// helper function to calculate stop residue of loop without overlap
core::Size
loop_start_without_overlap( core::pose::Pose const & pose, core::Size startres, core::Size const overlap );

core::kinematics::FoldTree
remove_all_jump_atoms( core::kinematics::FoldTree const & orig );

core::kinematics::FoldTree
remove_missing_jump_atoms( core::pose::Pose const & pose, core::kinematics::FoldTree const & orig );

/// @brief given a residue, rebuilds all missing atoms
void rebuild_missing_atoms( core::pose::Pose & pose, core::Size const resi );

/// @brief helper function that looks for the given residue in a fold tree and returns the jump that controls its 6D-DoFs
int find_jump_rec(
	core::kinematics::FoldTree const & ft,
	core::Size const residue );

/// @brief inserts the peptide edges to accomodate the new jump edge given
void insert_peptide_edges( core::kinematics::FoldTree & ft, core::kinematics::Edge const & jedge );

// @brief parses a string containing single integers and ranges. Returns a vector of all possible values.
utility::vector1< core::Size >
parse_length_string( std::string const & len_str );

/// @brief given a number 0 <= x < 1, calculate an integer M <= x <= N
/// NOTE THAT THIS FUNCTION MODIFIES THE PARAMETER
core::Size
extract_int( core::Real & num, core::Size const m, core::Size const n );

/// @brief counts the beta-bulges in the given segment. This simply counts all abego A's in the segment
core::Size
count_bulges( components::StructureData const & perm, std::string const & segment );

/// @brief dumps a pose into another pose as a new chain
void
add_chain_from_pose( core::pose::PoseCOP to_add, core::pose::PoseOP combined );

/// @brief adds residues from template_pose to pose.  If new_chain == true, creates covalent bond
void
add_residues_to_pose(
	core::pose::Pose & pose,
	core::pose::Pose const & template_pose,
	bool const new_chain );

core::Size
get_resid(
	components::StructureData const & sd,
	std::string const & resid_str );

/// @brief evaluate linear chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.  If cutpoint variants are present at chainbreak, will use
///  existing variants and not modify them.  If cutpoint variants are not
///  found will add them and then remove them once calculation is finished.
///  Eventually this should be merged into protocols/forge/methods/linear_chainbreak.cc
///  However, the changes I needed to make to that file break certain parts
///  of remodel
core::Real
linear_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos );

core::kinematics::FoldTree
slide_jump(
	core::kinematics::FoldTree const & ft_orig,
	core::Size const jump_idx,
	core::Size const new_start,
	core::Size const new_stop );

void
add_cutpoints( core::pose::Pose & pose, components::StructureData const & sd );

/// @brief Given a symmetric pose, and a fold tree for the asymmetric unit, constructs and
///        returns a symmetric fold tree while preserving the topology of the aysmmetric
///        unit's fold tree
core::kinematics::FoldTree
symmetric_fold_tree( core::pose::Pose const & pose, core::kinematics::FoldTree const & asymm_ft );

/// @brief Given a symmetric pose and a secstruct for the asymmetric unit, constructs and
///        returns a secondary structure string compatible with the symmetric pose based on
///        the secondary structure of the asymmetric unit
std::string
symmetric_secstruct( core::pose::Pose const & pose, std::string const & asymm_secstruct );

/// @brief Given a symmetric pose and a ResidueSubset for the asymmetric unit, constructs
///        and returns a residue subset compatible with the symmetric pose based on the
///        given asymmetric unit residue subset
core::select::residue_selector::ResidueSubset
symmetric_residue_subset( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset );

/// @brief Computes secondary structure string from the given motifs
/// @param[in]  motif_str  Motif string to be parsed (e.g. "5EB-2LG-5EB")
/// @param[out] secstruct  Secondary structure string to be cleared and filled
/// @param[out] abego      ABEGO string to be cleared and filled
void
parse_motif_string( std::string const & motif_str, std::string & secstruct, std::string & abego );

/// @brief Returns vector of elements found in v1 but not in v2
template< typename Container, typename T >
Container
set_difference(
	typename std::set< T >::const_iterator begin1,
	typename std::set< T >::const_iterator const & end1,
	typename std::set< T >::const_iterator begin2,
	typename std::set< T >::const_iterator const & end2
)
{
	Container result;
	std::set_difference( begin1, end1, begin2, end2, std::inserter( result, result.begin() ) );
	return result;
}

/// @brief Given a comma-seperated string, convert to container
template< typename Container >
Container
csv_to_container( std::string const & csv, char const delim=',' )
{
	utility::vector1< std::string > const fields = utility::string_split( csv, delim );
	return Container( fields.begin(), fields.end() );
}

} // denovo_design
} // protocols

#endif // INCLUDED_protocols_denovo_design_util_hh

