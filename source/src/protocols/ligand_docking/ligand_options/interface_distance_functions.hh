// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/ligand_docking/ligand_options/interface_distance_functions.hh
/// @brief Helper functions for the interface class
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_ligand_options_interface_distance_functions_hh
#define INCLUDED_protocols_ligand_docking_ligand_options_interface_distance_functions_hh

//// Utility Headers
#include <core/types.hh>

/// Project Headings

// STL Headers

#include <core/conformation/Residue.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {
namespace ligand_options {

bool check_all_ligand_atoms(
		core::conformation::Residue const & ligand_interface_residue,
		core::Vector const potential_interface_vector,
		core::Real const cutoff
);

bool check_neighbor_ligand_atom(
		core::conformation::Residue const & ligand_interface_residue,
		core::Vector const potential_interface_vector,
		core::Real const cutoff
);

bool is_interface_vector(
		core::Vector const & questionable_vector,
		core::Vector const & ligand_vector,
		core::Real const cutoff_squared
);

} //namespace ligand_options
} //namespace ligand_docking
} //namespace protocols

#endif
