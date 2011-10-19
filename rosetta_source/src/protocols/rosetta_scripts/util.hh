// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/RosettaScripts/util.hh
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_util_hh
#define INCLUDED_protocols_rosetta_scripts_util_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
// Utillity Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <set>

namespace protocols {
namespace rosetta_scripts {

core::Size
get_resnum( utility::tag::TagPtr const tag_ptr, core::pose::Pose const & pose, std::string const & prefix="" );

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
get_resnum_list(utility::tag::TagPtr const tag_ptr, std::string const& tag, core::pose::Pose const& pose);

/// @brief returns a resnum list directly from a string
std::set< core::Size >
get_resnum_list( std::string const str, core::pose::Pose const & pose );

/// @brief returns a resnum list directly from a string, preserving order
utility::vector1<core::Size>
get_resnum_list_ordered( std::string const str, core::pose::Pose const & pose );

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data );

core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data );

core::scoring::ScoreFunctionOP
parse_score_function( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data );


/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
core::pose::PoseOP
saved_reference_pose( utility::tag::TagPtr const in_tag, protocols::moves::DataMap & data_map );

void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm, protocols::moves::DataMap &, bool const reset_movemap = true /* should we turn everything to true at start?*/ );

void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm );

protocols::filters::FilterOP
parse_filter( std::string const filter_name, protocols::filters::Filters_map const & d );

protocols::moves::MoverOP
parse_mover( std::string const mover_name, protocols::moves::Movers_map const & d );

numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagPtr const xyz_vector_tag );

} // RosettaScripts
} // protocols


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
