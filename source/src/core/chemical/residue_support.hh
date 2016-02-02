// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/residue_support.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Steven Combs


#ifndef INCLUDED_core_chemical_residue_support_hh
#define INCLUDED_core_chemical_residue_support_hh

#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


namespace core {
namespace chemical {

// Find a better place to declare this function
/// @brief relies on class Graph to find all pairs shortest path information
ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res );

//used to create a light weight residue for searching rings
LightWeightResidueGraph convert_residuetype_to_light_graph(ResidueType const & res);

/// @brief Rename atoms in the residue type such that their names are unique.
/// If preserve is true, only rename those which have no names or who have
/// name conflicts. (Both of the conflicting atoms will be renamed.)
void
rename_atoms( ResidueType & res, bool preserve=true );


/// @brief Calculate the rigid matrix for neighbor atom finding
/// Assume that distances has been initialized to some really large value, and is square
void calculate_rigid_matrix( ResidueType const & res, utility::vector1< utility::vector1< core::Real > > & distances );

/// @brief Find the neighbor distance to the given neighbor atom.
/// If nbr_atom is null_vertex, give the smallest neighbor distance,
/// and set nbr_atom to the atom for that distance.
/// @details The neighbor distance here is adjusted for rotatable bonds -
/// It should be at least as large as the maximum neighbor distance
/// in any torsional rotamer
/// If the neighbor atom is not provided, the atom chosen will be a
/// multiply-bonded heavy atom.
///
/// Assumes:
///   * All atoms and bond are present
///   * All ideal_xyz coordinates have been set
///   * All elements have been set
///  * All ring bonds have been annotated
core::Real
find_nbr_dist( ResidueType const & res, VD & nbr_atom );

/// @brief Apply molfile_to_params style partial charges to the ResidueType.
/// @details These partial charges are based off of the Rosetta atom type,
/// adjusted such that the net partial charge is equal to the net formal charge.
///
/// These charges are almost certainly dodgy. If you have any other source of
/// partial charges that are at all reasonable, you probably want to consider those instead.
///
/// Assumes:
///   * All atoms and bond are present.
///   * All atom types have been set.
///   * Formal charges (if any) have been set.
void
rosetta_recharge_fullatom( ResidueType & res );


} // chemical
} // core

#endif
