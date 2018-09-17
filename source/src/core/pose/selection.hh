// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pose/selection.hh
/// @brief selcetion of pose residues
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu) willsheffler@gmail.com

#ifndef INCLUDED_core_pose_selection_hh
#define INCLUDED_core_pose_selection_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/ResidueIndexDescription.fwd.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace core {
namespace pose {

/// @brief Creates a ResidueIndexDescription from a string
/// @detail Recognizes three forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
///     IMPORTANT: This does not currently handle insertion codes.
///  - Reference pose numbers.  These have the form refpose([refpose name], [refpose number]).
/// In addition, relative numbers are permitted (of the form +[number] or -[number]) in conjunction
/// with reference pose numbering.  For example, one might say "refpose(state1,17)+3", which means
/// three residues past the residue correpsonding to residue 17 in the reference pose called "state1".
/// @return a ResidueIndexDescription that will yield a residue number when applied to pose.
///  Returns a nullptr if there's an error with parsing the string
ResidueIndexDescriptionCOP
parse_resnum(std::string const& resnum, bool const check_for_refpose=false);

/// @brief Extracts a residue number from a string.
/// @detail Recognizes three forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
///     IMPORTANT: This does not currently handle insertion codes.
///  - Reference pose numbers.  These have the form refpose([refpose name], [refpose number]).
/// In addition, relative numbers are permitted (of the form +[number] or -[number]) in conjunction
/// with reference pose numbering.  For example, one might say "refpose(state1,17)+3", which means
/// three residues past the residue correpsonding to residue 17 in the reference pose called "state1".
/// @return the rosetta residue number for the string, or 0 upon an error
core::Size
parse_resnum(std::string const& resnum, core::pose::Pose const& pose, bool const check_for_refpose=false);

/// @brief returns a resnum list directly from a string
std::set< core::Size >
get_resnum_list( std::string const & str, core::pose::Pose const & pose );

/// @brief returns a resnum list directly from a string, preserving order
utility::vector1<core::Size>
get_resnum_list_ordered( std::string const & str, core::pose::Pose const & pose );

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

/////////////////////////////////////////////////////////
//////////////////// XML PARSING ////////////////////////

/// @brief DEPRECATED - provided for legacy usage only. Don't use for new code.
/// Instead, just use a single option which uses the parse_resnum syntax to specify.
/// @details Checks for either/both of `pdb_num`/`res_num` in the tag, and pulls out the appropriate string
std::string
get_resnum_string( utility::tag::TagCOP tag_ptr, std::string const & prefix="" );

/// @brief DEPRECATED - provided for legacy usage only. Don't use for new code.
/// Instead, just use a single option which uses parse_resnum syntax to specify
/// @details Checks for either/both of `pdb_num`/`res_num` in the tag, and pulls out the appropriate string
std::string
get_resnum_string( utility::tag::TagCOP tag_ptr, std::string const & prefix, std::string const & default_value );

///// @brief Extracts a list of residue numbers from a tag
///// @details The tag should contain a comma-separated list of numbers, in either
/////   pdb or rosetta format (@see parse_resnum for details)
core::select::residue_selector::ResidueSelectorOP
get_resnum_selector(utility::tag::TagCOP tag_ptr, std::string const& tag);

//////////// XSD /////////////////

///@brief Companion function for get_resnum_string
///@details Appends relevant XMLSchemaAttributes to the AttributeList
void
attributes_for_get_resnum_string( utility::tag::AttributeList & attlist, std::string const & prefix="" );

///@brief Companion function for get_resnum_selector
///@details Appends relevant XMLSchemaAttributes to the AttributeList; needs the xsd because it adds restricted types
void
attributes_for_get_resnum_selector( utility::tag::AttributeList & attlist, utility::tag::XMLSchemaDefinition & xsd, std::string const& tag="" );

///@brief Companion function for parse_resnum
///@details Appends relevant XMLSchemaAttributes to the AttributeList
void
attributes_for_parse_resnum( utility::tag::AttributeList & attlist, std::string const & att_name, std::string const & description = "");

} // pose
} // core


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
