// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/ga_dock/util.hh
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_util_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_util_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <string>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

core::Size
count_neighbors_on_coord( core::pose::Pose const &pose,
	core::Vector const &xyz1,
	std::string const atomname,
	core::Real const dcut );

utility::vector1< core::Size >
count_neighbors( core::pose::Pose const &pose,
	std::string const atomname,
	core::Real const dcut );

utility::vector1< core::Size >
get_atomic_contacting_sidechains( core::pose::Pose const &pose,
	utility::vector1< core::Size > const &ligids,
	core::Real const atomic_distance_cut = 4.0 );

void
constraint_relax( core::pose::Pose &pose,
	utility::vector1< core::Size > const &ligids,
	utility::vector1< core::Size > const &movable_scs,
	core::Real maxiter = 50.0
);

void
make_ligand_only_pose( core::pose::PoseOP pose_new,
	core::pose::PoseCOP pose, //pass by value
	utility::vector1< core::Size > const& lig_resnos
);

void
perturb_ligand_rb(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const& ligids,
	core::Real trans_step = 2.0,
	core::Real rot_step = 15.0
);

void
perturb_ligand_torsions(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const& ligids,
	utility::vector1< core::Size > const& freeze_chi,
	core::Real chi_step = 15.0
);

core::Size
get_ligand_jumpid( core::pose::Pose const &pose,
	utility::vector1< core::Size > const& ligids );

}
}
}

#endif // INCLUDED_protocols_ligand_docking_util_HH
