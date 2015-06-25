// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/pose/selection.hh
/// @brief selcetion of pose residues
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu) willsheffler@gmail.com

#ifndef INCLUDED_core_pose_selection_hh
#define INCLUDED_core_pose_selection_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// Utillity Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace core {
namespace pose {

core::Size
get_resnum( utility::tag::TagCOP tag_ptr, core::pose::Pose const & pose, std::string const & prefix="" );

/// @brief Extracts a residue number from a string.
/// @detail Recognizes three forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
///		- Reference pose numbers.  These have the form refpose([refpose name], [refpose number]).
/// In addition, relative numbers are permitted (of the form +[number] or -[number]) in conjunction
/// with reference pose numbering.  For example, one might say "refpose(state1,17)+3", which means
/// three residues past the residue correpsonding to residue 17 in the reference pose called "state1".
/// @return the rosetta residue number for the string, or 0 upon an error
core::Size
parse_resnum(std::string const& resnum, core::pose::Pose const& pose, bool const check_for_refpose=false);

/// @brief Extracts a list of residue numbers from a tag.
/// @details The tag should contain a comma-separated list of numbers, in either
///   pdb or rosetta format (@see parse_resnum for details)
utility::vector1<core::Size>
get_resnum_list(utility::tag::TagCOP tag_ptr, std::string const& tag, core::pose::Pose const& pose);

/// @brief returns a resnum list directly from a string
std::set< core::Size >
get_resnum_list( std::string const &str, core::pose::Pose const & pose );

/// @brief returns a resnum list directly from a string, preserving order
utility::vector1<core::Size>
get_resnum_list_ordered( std::string const &str, core::pose::Pose const & pose );

/// @brief Is a string of the format "refpose(<refposename>,<refposenumber>)" or "refpose(<refposename>,<refposenumber>)+/-<number>"?
/// @details  If this successfully determines that this is a string of this format, it populates the refpose_string, refpose_resnumber,
/// and refpose_offset variables with the name of the ReferencePose, the number of the residue in the reference pose, and the +/- offset
/// number parsed from this string.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
bool is_referencepose_number(
	std::string const &str,
	std::string &refpose_string,
	core::Size &refpose_resnumber,
	signed long &refpose_offset
);

/// @brief Given the name of a ReferencePose object in the pose, a residue number in that reference pose, and a residue offset,
/// this function returns the Rosetta number of the corresponding residue in the pose.  Should throw an error if the ReferencePose
/// doesn't exist in the pose, or 0 if no corresponding residue exists in the pose.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
core::Size get_resnumber_from_reference_pose(
	std::string const &refpose_name,
	core::Size const refpose_number,
	signed long const refpose_offset,
	core::pose::Pose const &pose
);

} // pose
} // core


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
