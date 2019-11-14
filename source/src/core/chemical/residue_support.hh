// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/residue_support.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Steven Combs


#ifndef INCLUDED_core_chemical_residue_support_hh
#define INCLUDED_core_chemical_residue_support_hh

#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <core/chemical/Bond.fwd.hh>

namespace core {
namespace chemical {

// TODO: This function seems a bit out of place here
/// @brief relies on class Graph to find all pairs shortest path information
ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res );

/// @brief Figure out the shortest path between the upper and lower connect atoms (inclusive)
/// Will return an empty vector if one does not exist.
utility::vector1< VD >
mainchain_path( MutableResidueType const & res );

/// @brief Figure out the shortest path between two atoms (inclusive)
/// Will return an empty vector if one does not exist.
utility::vector1< VD >
shortest_path( MutableResidueType const & res, VD start, VD end );

/// @brief Annotate "backbone" atoms.
/// For the purpose of this function, backbone atoms are any atoms which are connected
/// to another backbone atom by a non-rotatable, non-cut bond.
/// Atoms connected to the upper and/or lower connect points are always backbone.
/// Important - if Chis/cuts aren't properly annotated, all atoms will be backbone.
void
annotate_backbone( MutableResidueType & restype );

/// @brief Virtualize convert the MutableResidueType to a virtual type
/// NOTE: This function does not rename the residue type
// RM: Not sure why this isn't a patch/patch operation
void
real_to_virtual( MutableResidueType & restype );

//used to create a light weight residue for searching rings
LightWeightResidueGraph convert_residuetype_to_light_graph(MutableResidueType const & res);

/// @brief Rename atoms in the residue type such that their names are unique.
/// If preserve is true, only rename those which have no names or who have
/// name conflicts. (Both of the conflicting atoms will be renamed.)
void
rename_atoms( MutableResidueType & res, bool preserve=true );

/// @brief Utility class for VD-indexed matrix
class VDDistanceMatrix {
private:
	typedef boost::unordered_map< VD, core::Real > InternalVector;
	typedef boost::unordered_map< VD, InternalVector > InternalMatrix;
public:

	VDDistanceMatrix(core::Real default_val):
		default_(default_val)
	{}

	core::Real &
	operator() ( VD a, VD b );

	core::Real
	find_max_over( VD a );

private:

	InternalMatrix matrix_;
	core::Real default_;

};

/// @brief Calculate the rigid matrix - assume that distances has been initialized to some really large value, and is square
/// Helper for find_nbr_dist
void calculate_rigid_matrix( MutableResidueType const & res, VDDistanceMatrix & distances );

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
find_nbr_dist( MutableResidueType const & res, VD & nbr_atom );

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
rosetta_recharge_fullatom( MutableResidueType & res );

/// @brief Make a centroid version of the fullatom ResidueType passed in.
///
/// May return a nullptr if the conversion is not possible
///
/// @details This uses the same crude heuristics as in molfile_to_params
/// That is, all heavy atoms and polar hydrogens are transfered over 1:1, and
/// non-polar hydrogens are deleted.
///
/// Of particular note is that it makes no attempt to transfer things over into
/// "Superatoms"
///
/// Current limitation: it cannot convert a ResidueType which has connections,
/// if any of the ICOORs depend on deleted hydrogens.
///
/// Assumes:
///   * Input ResidueType is complete and finalized
MutableResidueTypeOP
make_centroid( ResidueType const & res );

MutableResidueTypeOP
make_centroid( MutableResidueType const & res );

/// @brief Are two ResidueTypes equivalent?
/// This is here rather than as an operator on ResidueType because it's not the sort of thing one should be doing normally.
/// This looks for exact equivalence, including atom order.
bool
residue_types_identical( ResidueType const & res1, ResidueType const & res2 );

/// @brief Are the two ResidueConnection objects equivalent
/// Here instead of in operator== because of the fuzzy-real issue.
bool
compare_residue_connection( ResidueConnection const & rc1, ResidueConnection const & rc2, bool fuzzy=false );

/// @brief Are the two ResidueConnection objects equivalent
/// Here instead of in operator== because of the fuzzy-real issue.
bool
compare_atom_icoor( AtomICoor const & aic1, AtomICoor const & aic2, bool fuzzy=false );

} // chemical
} // core

#endif
