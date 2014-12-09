// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/util.hh
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_constraints_UTIL_HH
#define INCLUDED_protocols_antibody_constraints_UTIL_HH

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyEnum.hh>



namespace protocols {
namespace antibody {
namespace constraints {
	
///@brief Check if all residues already have a specific constraint type.  Useful for coordinate and dihedral constraints.
bool
cdr_has_res_constraints(AntibodyInfoCOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, std::string const constraint_type);


///@brief Adds harmonic cluster constraint or coordinate cst if the cluster does not have constraints associated with it. 
/// Requires AHO numbering.  
void
add_harmonic_cluster_cst_or_coordinate_cst(AntibodyInfoOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, core::Real coord_cst_sd = .75);
	
///@brief Adds harmonic cluster constraint or coordinate cst if the cluster does not have constraints associated with it. 
/// Requires AHO numbering.  
void
add_harmonic_cluster_cst_or_coordinate_cst(AntibodyInfoOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, clusters::CDRClusterEnum const cluster, core::Real coord_cst_sd = .75);


///@brief Add dihedral constraints to CDR with mean being the current phi/psi if no cluster constraints.
///@details Values for avg SD for each cluster not including H3 were 23 and 42 degrees respectively.  Values chosen here are conservative.
///                   Mean SD for both dihedrals was 32.66 degrees
///                   Use this especially if doing cartesian-space minimization
///
void
add_harmonic_cluster_cst_or_dihedral_cst(
	AntibodyInfoOP ab_info, core::pose::Pose & pose, 
	CDRNameEnum const cdr,
	core::Real phi_sd_deg = 20.0, core::Real psi_sd_deg = 30.0);

///@brief Add dihedral constraints to CDR with mean being the current phi/psi if no cluster constraints.
///@details Values for avg SD for each cluster not including H3 were 23 and 42 degrees respectively.  Values chosen here are conservative.
///                   Mean SD for both dihedrals was 32.66 degrees
///                   Use this especially if doing cartesian-space minimization
///
void
add_harmonic_cluster_cst_or_dihedral_cst(
	AntibodyInfoOP ab_info, core::pose::Pose & pose, 
	CDRNameEnum const cdr, clusters::CDRClusterEnum const cluster,
	core::Real phi_sd_deg = 20.0, core::Real psi_sd_deg = 30.0);


///@brief Add dihedral constraints to CDR with mean being the current phi/psi.
///@details Values for avg SD for each cluster not including H3 were 23 and 42 degrees respectively.  Values chosen here are conservative.
///                   Mean SD for both dihedrals was 32.66 degrees
///                   Use this especially if doing cartesian-space minimization
///
void
add_harmonic_dihedral_cst_general(
	AntibodyInfoOP ab_info, core::pose::Pose & pose,
	CDRNameEnum const cdr,
	core::Real phi_sd_deg = 20.0, core::Real psi_sd_deg = 30.0);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North CDR Cluster Constraints. 
//                     Requires AHO numbering due to hardcoded constraint files.
//
//


/// @brief Adds dihedral harmonic constraints to Pose CDRs using cluster info in AntibodyInfo
/// @details Currently requires North_AHO numbering. Returns map of success/failure
std::map<CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, core::pose::Pose & pose);

///@brief Same as above, but adds constraints to the vector so they can be identified and removed from the pose if needed.
///@details Returns map of success/failure
std::map<CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, core::pose::Pose & pose, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints);


/// @brief Adds a harmonic constraint to a Pose CDR based on cluster type
/// @details Currently requires North_AHO numbering.
bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, core::pose::Pose & pose, clusters::CDRClusterEnum const cluster);

///@brief Same as above, but adds constraints to the vector so they can be identified and removed from the pose if needed.
///@details Returns true or false depending on success
bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, core::pose::Pose & pose, clusters::CDRClusterEnum const cluster, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints);

/// @brief Gets the cluster constraint name.  Returns NA if not found.
std::string
get_harmonic_cluster_constraint_filename(AntibodyInfoCOP ab_info, clusters::CDRClusterEnum const cluster);


}
}
}

#endif	//#ifndef INCLUDED_protocols/antibody_design_UTIL_HH
