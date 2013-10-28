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
get_resnum( utility::tag::TagCOP const tag_ptr, core::pose::Pose const & pose, std::string const & prefix="" );

/// @brief Extracts a residue number from a string.
/// @detail Recognizes two forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
/// @return the rosetta residue number for the string, or 0 upon an error
core::Size
parse_resnum(std::string const& resnum, core::pose::Pose const& pose);

/// @brief Extracts a list of residue numbers from a tag.
/// @details The tag should contain a comma-separated list of numbers, in either
///   pdb or rosetta format (@see parse_resnum for details)
utility::vector1<core::Size>
get_resnum_list(utility::tag::TagCOP const tag_ptr, std::string const& tag, core::pose::Pose const& pose);

/// @brief returns a resnum list directly from a string
std::set< core::Size >
get_resnum_list( std::string const str, core::pose::Pose const & pose );

/// @brief returns a resnum list directly from a string, preserving order
utility::vector1<core::Size>
get_resnum_list_ordered( std::string const str, core::pose::Pose const & pose );


} // pose
} // core


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
