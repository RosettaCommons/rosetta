// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/jump_from_internal_coords.hh
/// @brief  Utility function for calculating jumps by knowing desired internal coordinates of arbitrary atoms
/// @author Jack Maguire


#ifndef INCLUDED_core_kinematics_inverse_jump_from_internal_coords_HH
#define INCLUDED_core_kinematics_inverse_jump_from_internal_coords_HH


// Package headers
#include <core/id/types.hh>

#include <core/types.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Conformation.hh>

#include <numeric/xyzVector.hh>

#include <array>

namespace core {
namespace kinematics {
namespace inverse {

struct InternalCoordAtoms {
	id::AtomID child;
	id::AtomID parent;
	id::AtomID grandparent;
};

struct InternalCoordGeometry {
	// NOTATION:

	// grandparent
	// | parent
	// | | child
	// | | |
	// | | |   child
	// | | |   | parent
	// | | |   | | grandparent
	// | | |   | | |
	// A-B-C---D-E-F
	// ^---^   ^---^
	// Fixed   Moving

	//calculate D's position
	Real A_B_C_D_torsion_angle_rad = 0; //radians
	Real B_C_D_bond_angle_rad = 0; //radians - IMPROPER BOND ANGLE
	Real C_D_dist_Ang = 0; //Angstroms

	//calculate E's position
	Real B_C_D_E_torsion_angle_rad = 0; //radians
	Real C_D_E_bond_angle_rad = 0; //radians - IMPROPER BOND ANGLE

	//calculate F's position
	Real C_D_E_F_torsion_angle_rad = 0; //radians


	///@brief set all six values using their current configuration in the conformation
	void
	init_from_current(
		conformation::Conformation const & conformation,
		InternalCoordAtoms const & ABC,
		InternalCoordAtoms const & DEF
	);
};

///@brief Utility function for calculating jumps by knowing desired internal coordinates of arbitrary atoms
///@details fixed_atoms (A/B/C) are upstream of the jump, moving_atoms (D/E/F) are downstream.
///The user provides distances, angles, and dihedral angles for the relationship between the two sets and
///the number of the jump the user wants to change. Rosetta then calculates the new value of this jump.
Jump
jump_from_internal_coords(
	conformation::Conformation const & conformation,
	InternalCoordAtoms const & fixed_atoms,
	InternalCoordAtoms const & moving_atoms,
	InternalCoordGeometry geometry,
	Size jump_id
);

///@brief Given internal coordinates and XYZs of the parents, calculate XYZ of a child
numeric::xyzVector< Real >
calc_new_atom_location(
	numeric::xyzVector< Real > const & greatgrandparent,
	numeric::xyzVector< Real > const & grandparent,
	numeric::xyzVector< Real > const & parent,
	Real dist,         //angstroms
	Real bond_angle,   //radians
	Real torsion_angle //radians
);


#ifdef PYROSETTA
// Force creation of PyRosetta bindings for internal type(s)

/*
From Sergey in PR #5634:
The easiest way to trigger binding generation for this type would be to:
-- in header file that uses this type
-- create #ifdef PYROSETTA block
-- define inline function that that this type by-value
*/
inline
InternalCoordAtoms
dummy_InternalCoordAtoms( InternalCoordAtoms val ){
	return val;
}

inline
InternalCoordGeometry
dummy_InternalCoordGeometry( InternalCoordGeometry val ){
	return val;
}
#endif

} // namespace inverse
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_inverse_jump_from_internal_coords_HH
