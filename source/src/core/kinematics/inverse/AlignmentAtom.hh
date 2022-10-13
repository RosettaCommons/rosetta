// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/AlignmentAtom.hh
/// @brief  Utility functions for calculating jumps by knowing desired atom positions
/// @author Jack Maguire


#ifndef INCLUDED_core_kinematics_inverse_AlignmentAtom_HH
#define INCLUDED_core_kinematics_inverse_AlignmentAtom_HH

// Package headers
#include <core/kinematics/inverse/AlignmentAtom.fwd.hh>

#include <core/id/types.hh>

#include <core/types.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/conformation/Conformation.hh>

#include <numeric/xyzVector.hh>

#include <array>

namespace core {
namespace kinematics {
namespace inverse {

///@brief this class has the duty of logging atoms and their desired positions
///@details We are using struct instead of class because AlignmentAtom has no invariants
struct AlignmentAtom {
	AlignmentAtom(
		core::id::AtomID atomid,
		numeric::xyzVector< core::Real > const dest_xyz
	):
		id( atomid ),
		destination_xyz( dest_xyz )
	{}

	AlignmentAtom() = default;

	void set(
		core::id::AtomID atomid,
		numeric::xyzVector< core::Real > const dest_xyz
	){
		id = atomid;
		destination_xyz = dest_xyz;
	}

	void set(
		core::id::AtomID atomid,
		conformation::Conformation const & conf
	){
		set( atomid, conf.xyz( atomid ) );
	}

	core::id::AtomID id;
	numeric::xyzVector< core::Real > destination_xyz;
};

///@brief This class contains the 3 AlignmentAtoms that we want to use to define the new jump value.
///@details In theory, the order of the atoms does not matter. We are using struct instead of class because AlignmentAtomArray has no invariants. We are also defining functions in the header because the function body does a better job of explaining the purpose than my comments do.
struct AlignmentAtomArray {

	std::array< AlignmentAtom, 3 > atoms;

	///@brief Define the intended stub that represents the future pose _after_ updating the jump
	Stub
	create_destination_stub() const {
		return Stub(
			atoms[0].destination_xyz,
			atoms[1].destination_xyz,
			atoms[2].destination_xyz
		);
	}

	///@brief This utility function uses the current XYZ coordinates for the AtomIDs in the pose.
	Stub
	create_stub_from_atom_ids( conformation::Conformation const & conf ) const {
		return Stub(
			conf.xyz( atoms[0].id ),
			conf.xyz( atoms[1].id ),
			conf.xyz( atoms[2].id )
		);
	}

}; //struct AlignmentAtomArray


#ifdef PYROSETTA
// Force creation of PyRosetta bindings for internal type(s)

//This is just future-proofing in case someone changes the type above but doesn't update down here
using AlignmentAtomArrayInternalType = decltype( AlignmentAtomArray::atoms );

/*
From Sergey in PR #5634:
The easiest way to trigger binding generation for this type would be to:
-- in header file that uses this type
-- create #ifdef PYROSETTA block
-- define inline function that that this type by-value
*/
inline
AlignmentAtomArrayInternalType
dummy_AlignmentAtomArrayInternalType( AlignmentAtomArrayInternalType val ){
	return val;
}
#endif

} // namespace inverse
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_inverse_AlignmentAtom_HH
