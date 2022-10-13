// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

///@file   core/kinematics/inverse/jump.cc
///@brief  Utility functions for calculating jumps by knowing desired atom positions
///@author Jack Maguire

// Unit headers
#include <core/kinematics/inverse/jump.hh>
#include <core/kinematics/inverse/util.hh>
#include <core/kinematics/inverse/AlignmentAtom.hh>

#include <basic/Tracer.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <utility/excn/Exceptions.hh>
//#include <stdexcept>

namespace core {
namespace kinematics {
namespace inverse {

static basic::Tracer TR( "core.kinematics.inverse.jump" );

using XYZ = numeric::xyzVector< core::Real >;

///@brief Update a desired jump in the conformation to place the atoms in the AlignmentAtomArray to the desired position
///@details protocol is described as comments in the body of the code.
///@author Jack Maguire
Jump
calculate_new_jump(
	conformation::Conformation const & conformation,
	core::Size const jump_id,
	AlignmentAtomArray const & atoms
){
	assert_atoms_are_downstream_of_jump( conformation, jump_id, atoms );
	// Notation:
	// Stub "A" is the root of the jump, which goes unchanged
	// Stub "B" is the destination of the jump, which we're hoping to change. B0 is the input position, B1 is the new position
	// Stub "C" is made of the 3 atoms being used to define our desired outcome. C0 is the starting position, C1 is the desired position as defined by the user

	// jump_A_B connects Stub A to Stub B (jump_A_B0 goes to stub B0)
	// Note, jump_B0_C0 is equal to jump_B1_C1 so we're just going to call that constant jump_B_C
	// jump_A_B will have an anologous homogeneous transform ht_A_B

	// Current Situation:
	// stubA -> stubB0 -> stubC0
	//       ^         ^
	//     jump_A_B0   jump_B_C

	// Desired Situation
	// stubA -> stubB1 -> stubC1
	//       ^         ^
	//     jump_A_B1   jump_B_C


	// We want to solve for ht_A_B1 (which maps to jump_A_B1) using the equation:
	// ht_A_B1 * ht_B_C = ht_A_C1
	// ht_A_B1 = ht_A_C1 / ht_B_C


	/////////////
	//kinematics

	//stubs
	Stub const & stub_A = conformation.upstream_jump_stub( jump_id );
	Stub const & stub_B0 = conformation.downstream_jump_stub( jump_id );
	Stub const & stub_C0 = atoms.create_stub_from_atom_ids( conformation );
	Stub const & stub_C1 = atoms.create_destination_stub();

	//jumps
	Jump const jump_B_C( stub_B0, stub_C0 );
	Jump const jump_A_C1( stub_A, stub_C1 );

	//homogeneous transforms
	using HT = numeric::HomogeneousTransform< core::Real >;
	HT const ht_B_C(   jump_B_C.get_rotation(),  jump_B_C.get_translation() );
	HT const ht_A_C1( jump_A_C1.get_rotation(), jump_A_C1.get_translation() );
	HT const ht_A_B1 = ht_A_C1 * ht_B_C.inverse();

	Jump jump_A_B1 = conformation.jump( jump_id );
	jump_A_B1.set_rotation( ht_A_B1.rotation_matrix() );
	jump_A_B1.set_translation( ht_A_B1.point() );

	return jump_A_B1;

	//Can do:
	//pose.set_jump( jump_id, jump_A_B1 );
}



} // namespace inverse
} // namespace kinematics
} // namespace core
