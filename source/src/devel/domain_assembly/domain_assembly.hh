// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_devel_domain_assembly_domain_assembly_hh
#define INCLUDED_devel_domain_assembly_domain_assembly_hh

#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/farna/RNA_FragmentsClasses.fwd.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

namespace devel {
namespace domain_assembly {


/// @brief optimizes linkers in a multidomain protein
void assemble_domains_optimize();

/// @brief Initialize a movemap reading from the command line option
/// -da_linker_file, loading linker residue ranges from this file, and
/// then initializing the movemap from these ranges
core::kinematics::MoveMapOP
read_movemap_from_da_linker_file();

/// @brief reads in file that specifies which regions of the protein will
///  move during domain assembly
///  Each line of the file should have the start and end position for a linker region
bool
read_linker_file(
	std::string const & filename,
	utility::vector1< std::pair < core::Size, core::Size > > & linker_ranges
);

/// @brief sets movemap true for regions specified in linker file
void
set_movemap_for_linkers(
	utility::vector1< std::pair < core::Size, core::Size > > const & linker_ranges,
	core::kinematics::MoveMapOP & mm
);

/// @brief centroid mode optimization of linkers
void
optimize_linkers_centroid_mode(
	core::kinematics::MoveMapOP & mm,
	core::pose::Pose & full_pose
);

void
optimize_linkers_fullatom_mode(
	core::kinematics::MoveMapOP & mm,
	core::pose::Pose & full_pose
);

ObjexxFCL::FArray1D_bool
set_moveable_rna(
	core::pose::Pose & full_pose,
	utility::vector1< std::pair < core::Size, core::Size > > & linker_rna
);

void
optimize_linkers_rna_fullatom_mode(
	core::kinematics::MoveMapOP & mm,
	core::pose::Pose & full_pose,
	protocols::farna::RNA_FragmentsOP & all_rna_fragments
);

/// @brief a helper function for the domain assembly protocol. Selects
///residues near linkers and domain interfaces for repacking
void
da_residues_to_repack(
	core:: kinematics::MoveMapOP & mm,
	utility::vector1< std::pair < core::Size, core::Size > > & nearest_movable_residues,
	core::pose::Pose & pose,
	utility::vector1<bool> & repack_residues
);

/// @brief a helper function for the domain assembly protocol.  For each residue
/// it finds the closest N-terminal and C-terminal movable residue (as specified
/// in the input movemap)
void
find_nearest_movable_residues(
	core::kinematics::MoveMapOP & mm,
	core::pose::Pose & pose,
	utility::vector1< std::pair < core::Size, core::Size > > & nearest_movable_residues
);


}
}

#endif
