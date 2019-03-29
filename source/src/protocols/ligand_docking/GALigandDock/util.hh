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

#include <core/pose/Pose.hh>
#include <core/types.hh>

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
	core::Size const ligid,
	core::Real const atomic_distance_cut = 4.0 );

}
}
}

#endif // INCLUDED_protocols_ligand_docking_util_HH
