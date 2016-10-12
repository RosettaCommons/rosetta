// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/util.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_UTIL_HH
#define INCLUDED_protocols_antibody_clusters_UTIL_HH

#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace clusters {

void
add_cluster_comments_to_pose(core::pose::Pose& pose, AntibodyInfoCOP ab_info);

std::string
get_pose_cis_trans_conformation(core::pose::Pose const & pose, core::Size const start, core::Size const end);

/// @brief Calculates the dihedral distance used to match cluster centers.
core::Real
calculate_dihedral_distance(
	utility::vector1< core::Real > cluster_phis,
	utility::vector1< core::Real > pose_phis,
	utility::vector1< core::Real > cluster_psis,
	utility::vector1< core::Real > pose_psis
);

/// @brief Very basic way to check to make sure pose residues are North_AHO (North, et al) scheme, which allows the clustering.
/// @details If any of these anchor residues that are checked are missing, it will return false.
bool
check_if_pose_renumbered_for_clusters(core::pose::Pose const & pose);

CDRClusterEnum
get_cluster_from_cache_or_ab_info(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, CDRNameEnum const cdr);

} //clusters
} //antibody
} //protocols


#endif //#ifndef INCLUDED_protocols/antibody/clusters/UTIL_HH

