// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/util.hh
/// @brief Utility functions for the recon_design namespace
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_recon_design_util_hh
#define INCLUDED_protocols_recon_design_util_hh

#ifdef USEMPI
#include <mpi.h>
#endif

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace recon_design {

/// @brief Based on a pose and a resfile, get the indices of all designable residues
utility::vector1< core::Size >
get_designable_residues( core::pose::Pose & pose, std::string resfile );

/// @brief Based on a pose and the indices of all designable residues, get a
/// string of all designable AAs concatenated
utility::vector1< std::string >
get_designable_sequence ( core::pose::Pose & pose, utility::vector1< core::Size > designable_residues );

/// @brief Based on a list of sequences from poses, get all the AAs present at
/// the position given by position_no
utility::vector1< std::string > get_candidate_AAs(
	utility::vector1< utility::vector1< std::string > > const & other_pose_sequences,
	core::Size position_no );

/// @brief Given a list of poses, find the index of a particular pose
core::Size find_pose_in_vector( core::pose::Pose const & pose, utility::vector1< core::pose::PoseOP > & poses );

} //namespace recon_design
} //namespace protocols

#endif //INCLUDED_protocols_recon_design_util_HH




