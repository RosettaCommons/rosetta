// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/ligand_docking/ligand_options/interface_distance_functions.cc
/// @brief  helper functions for the Interface class
/// @author Gordon Lemmon (glemmon@gmail.com)
///
// Unit Headers
#include <protocols/ligand_docking/ligand_options/interface_distance_functions.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {
namespace ligand_options {

static THREAD_LOCAL basic::Tracer interface_distance_tracer( "protocols.ligand_docking.ligand_options.InterfaceBuilder", basic::t_debug );

bool check_all_ligand_atoms(
	core::conformation::Residue const & ligand_interface_residue,
	core::Vector const & potential_interface_vector,
	core::Real const cutoff
){
	core::Real cutoff_squared= cutoff*cutoff;
	for ( core::Size i = 1, i_end = ligand_interface_residue.nheavyatoms(); i <= i_end; ++i ) {
		core::Vector ligand_vector= ligand_interface_residue.xyz(i);
		if ( is_interface_vector(potential_interface_vector, ligand_vector, cutoff_squared) ) {
			return true;
		}
	}
	return false;
}

bool check_neighbor_ligand_atom(
	core::conformation::Residue const & ligand_interface_residue,
	core::Vector const & potential_interface_vector,
	core::Real const cutoff
){
	core::Real cutoff_squared= cutoff*cutoff;
	core::Size ligand_neighbor_atom_id= ligand_interface_residue.nbr_atom();
	core::Vector ligand_vector= ligand_interface_residue.xyz(ligand_neighbor_atom_id);

	return is_interface_vector(potential_interface_vector, ligand_vector, cutoff_squared);
}

bool is_interface_vector(
	core::Vector const & questionable_vector,
	core::Vector const & ligand_vector,
	core::Real const cutoff_squared
){

	double distance_squared = ligand_vector.distance_squared(questionable_vector);
	if ( distance_squared <= cutoff_squared ) {
		return true; // don't need to go through all atoms
	} else {
		return false;
	}
}

} //namespace ligand_options
} //namespace ligand_docking
} //namespace protocols
