// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/util.hh
/// @brief
/// @author Colin A. Smith (colin.smith@ucsf.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) (copied from app)

#ifndef INCLUDED_protocols_backrub_util_hh
#define INCLUDED_protocols_backrub_util_hh


#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace backrub {

bool
read_fold_tree_from_file( core::kinematics::FoldTree & foldtree, std::string filepath );

bool
read_fold_tree_from_file( core::pose::Pose & pose, std::string filepath );

void
append_fold_tree_to_file( core::kinematics::FoldTree const & foldtree, std::string file_path );

utility::vector1<core::Size>
positions_incompatible_with_task( core::pose::Pose & pose, core::pack::task::PackerTask & packertask );

utility::vector1<core::Size>
get_pivot_residues_from_movemap( core::kinematics::MoveMapCOP movemap);


} //backrub
} //protocols

#endif //#ifndef INCLUDED_protocols_backrub_UTIL_HH

