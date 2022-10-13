// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

///@file   core/kinematics/inverse/jump_from_internal_coords.cc
///@brief  Utility functions for calculating jumps by knowing desired jump_from_internal_coords measurements
///@author Jack Maguire

// Unit headers
#include <core/kinematics/inverse/jump_from_internal_coords.hh>
#include <core/kinematics/inverse/jump.hh>
#include <core/kinematics/inverse/util.hh>
#include <core/kinematics/inverse/AlignmentAtom.hh>

#include <basic/Tracer.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/excn/Exceptions.hh>

namespace core {
namespace kinematics {
namespace inverse {

static basic::Tracer TR( "core.kinematics.inverse.jump_from_internal_coords" );

using XYZ = numeric::xyzVector< Real >;

void
InternalCoordGeometry::init_from_current(
	conformation::Conformation const & conf,
	InternalCoordAtoms const & ABC,
	InternalCoordAtoms const & DEF
){
	id::AtomID const A = ABC.grandparent;
	id::AtomID const B = ABC.parent;
	id::AtomID const C = ABC.child;
	id::AtomID const D = DEF.child;
	id::AtomID const E = DEF.parent;
	id::AtomID const F = DEF.grandparent;

	XYZ const & A_xyz = conf.xyz( A );
	XYZ const & B_xyz = conf.xyz( B );
	XYZ const & C_xyz = conf.xyz( C );
	XYZ const & D_xyz = conf.xyz( D );
	XYZ const & E_xyz = conf.xyz( E );
	XYZ const & F_xyz = conf.xyz( F );

	A_B_C_D_torsion_angle_rad = numeric::dihedral_radians( A_xyz, B_xyz, C_xyz, D_xyz );
	B_C_D_bond_angle_rad = numeric::angle_of( C_xyz-D_xyz, C_xyz-B_xyz );
	C_D_dist_Ang = C_xyz.distance( D_xyz );

	B_C_D_E_torsion_angle_rad = numeric::dihedral_radians( B_xyz, C_xyz, D_xyz, E_xyz );
	C_D_E_bond_angle_rad = numeric::angle_of( D_xyz-E_xyz, D_xyz-C_xyz );

	C_D_E_F_torsion_angle_rad = numeric::dihedral_radians( C_xyz, D_xyz, E_xyz, F_xyz );
}


///@brief Utility function for calculating jumps by knowing desired internal coordinates of arbitrary atoms
///@details fixed_atoms (A/B/C) are upstream of the jump, moving_atoms (D/E/F) are downstream.
///The user provides distances, angles, and dihedral angles for the relationship between the two sets and
///the number of the jump the user wants to change. Rosetta then calculates the new value of this jump.
Jump
jump_from_internal_coords(
	conformation::Conformation const & conf,
	InternalCoordAtoms const & fixed_atoms,
	InternalCoordAtoms const & moving_atoms,
	InternalCoordGeometry geom, //intentional copy
	Size const jump_id
) {
	runtime_assert( geom.C_D_dist_Ang > 0 );

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

	id::AtomID const A = fixed_atoms.grandparent;
	id::AtomID const B = fixed_atoms.parent;
	id::AtomID const C = fixed_atoms.child;
	id::AtomID const D = moving_atoms.child;
	id::AtomID const E = moving_atoms.parent;
	id::AtomID const F = moving_atoms.grandparent;

	XYZ const & A_xyz = conf.xyz( A );
	XYZ const & B_xyz = conf.xyz( B );
	XYZ const & C_xyz = conf.xyz( C );
	XYZ const & D_xyz = conf.xyz( D );
	XYZ const & E_xyz = conf.xyz( E );
	XYZ const & F_xyz = conf.xyz( F );

	///////////////
	//SANITY CHECKS
	{
		AlignmentAtomArray aaa;
		aaa.atoms[ 0 ].set( A, A_xyz );
		aaa.atoms[ 1 ].set( B, B_xyz );
		aaa.atoms[ 2 ].set( C, C_xyz );
		assert_atoms_are_upstream_of_jump( conf, jump_id, aaa );
	}
	{
		AlignmentAtomArray aaa;
		aaa.atoms[ 0 ].set( D, D_xyz );
		aaa.atoms[ 1 ].set( E, E_xyz );
		aaa.atoms[ 2 ].set( F, F_xyz );
		assert_atoms_are_downstream_of_jump( conf, jump_id, aaa );
	}


	///////////////////////////
	//ADDITIONAL GEOMETRY CALCS

	//Our bond angles are actually IMPROPER bond angles and need to be subtracted from pi
	using numeric::constants::d::pi;
	geom.B_C_D_bond_angle_rad = pi - geom.B_C_D_bond_angle_rad;
	geom.C_D_E_bond_angle_rad = pi - geom.C_D_E_bond_angle_rad;

	//Make angles fall within [0, 2*pi]
	while ( geom.B_C_D_bond_angle_rad < 0.0  ) geom.B_C_D_bond_angle_rad += 2*pi;
	while ( geom.B_C_D_bond_angle_rad > 2*pi ) geom.B_C_D_bond_angle_rad -= 2*pi;
	if ( geom.B_C_D_bond_angle_rad > pi ) {
		// larger-than-pi angles play funny.
		//We need to flip the corresponding torsion angle if we fall into this range

		// I can't fully explain why we don't need to change geom.B_C_D_bond_angle_rad itself but
		//there is a unit test to ensure that the expected behavior occurs
		geom.B_C_D_E_torsion_angle_rad += pi;
	}

	while ( geom.C_D_E_bond_angle_rad < 0.0  ) geom.C_D_E_bond_angle_rad += 2*pi;
	while ( geom.C_D_E_bond_angle_rad > 2*pi ) geom.C_D_E_bond_angle_rad -= 2*pi;
	if ( geom.C_D_E_bond_angle_rad > pi ) {
		geom.C_D_E_F_torsion_angle_rad += pi;
	}

	//Keep current geometries for unchanging angles
	Real const D_E_dist_Ang = D_xyz.distance( E_xyz );
	Real const E_F_dist_Ang = E_xyz.distance( F_xyz );
	Real const D_E_F_bond_angle_rad = numeric::angle_of( E_xyz-D_xyz, E_xyz-F_xyz );

	/////////////////////////////
	//CALCULATE NEW XYZ POSITIONS
	XYZ const new_D_xyz = calc_new_atom_location(
		A_xyz,
		B_xyz,
		C_xyz,
		geom.C_D_dist_Ang,
		geom.B_C_D_bond_angle_rad,
		geom.A_B_C_D_torsion_angle_rad
	);

	XYZ const new_E_xyz = calc_new_atom_location(
		B_xyz,
		C_xyz,
		new_D_xyz,
		D_E_dist_Ang,
		geom.C_D_E_bond_angle_rad,
		geom.B_C_D_E_torsion_angle_rad
	);

	XYZ const new_F_xyz = calc_new_atom_location(
		C_xyz,
		new_D_xyz,
		new_E_xyz,
		E_F_dist_Ang,
		D_E_F_bond_angle_rad,
		geom.C_D_E_F_torsion_angle_rad
	);

	////////////////////
	//CALCULATE NEW JUMP
	AlignmentAtomArray aaa;
	aaa.atoms[ 0 ].set( D, new_D_xyz );
	aaa.atoms[ 1 ].set( E, new_E_xyz );
	aaa.atoms[ 2 ].set( F, new_F_xyz );

	//offload the logic to kinematics/inverse/jump.hh
	return calculate_new_jump( conf, jump_id, aaa );
}

///@brief Given internal coordinates and XYZs of the parents, calculate XYZ of a child
numeric::xyzVector< Real >
calc_new_atom_location(
	numeric::xyzVector< Real > const & greatgrandparent,
	numeric::xyzVector< Real > const & grandparent,
	numeric::xyzVector< Real > const & parent,
	Real const dist,         //angstroms
	Real const bond_angle,   //radians
	Real const torsion_angle //radians
) {
	using namespace core::kinematics;
	using namespace core::kinematics::tree;

	std::array< tree::BondedAtomOP, 4 > atoms;

	for ( Size i = 0; i < 4; ++i ) {
		atoms[ i ] = utility::pointer::make_shared< BondedAtom >();
	}

	atoms[ 0 ]->append_atom( atoms[ 1 ] );
	atoms[ 1 ]->append_atom( atoms[ 2 ] );
	atoms[ 2 ]->append_atom( atoms[ 3 ] );

	atoms[ 0 ]->xyz( greatgrandparent );
	atoms[ 1 ]->xyz( grandparent );
	atoms[ 2 ]->xyz( parent );

	atoms[ 0 ]->update_internal_coords( true );

	atoms[ 3 ]->set_dof( id::PHI, torsion_angle );
	atoms[ 3 ]->set_dof( id::THETA, bond_angle );
	atoms[ 3 ]->set_dof( id::D, dist );

	atoms[ 3 ]->update_xyz_coords();

	return atoms[ 3 ]->xyz();
}

} // namespace inverse
} // namespace kinematics
} // namespace core

