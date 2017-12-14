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


/// @brief Check if all residues already have a specific constraint type.  Useful for coordinate and dihedral constraints.
bool
cdr_has_res_constraints(AntibodyInfoCOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, std::string const & constraint_type);


/// @brief Add dihedral constraints to CDR with mean being the current phi/psi.
/// @details Values for avg SD for each cluster not including H3 were 16 and 16 degrees respectively.  Using cluster outliers, values were at 21 and 23 degrees respectively.
///                   Mean SD for both dihedrals was 16 and 22 for without/with outliers.
///                   Use this especially if doing cartesian-space minimization
///
void
add_harmonic_dihedral_cst_general(
	AntibodyInfoCOP ab_info, core::pose::Pose & pose,
	CDRNameEnum const cdr,
	core::Real phi_sd_deg = 16.0, core::Real psi_sd_deg = 16.0);

}
}
}

#endif //#ifndef INCLUDED_protocols/antibody_design_UTIL_HH
