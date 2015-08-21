// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/backrub/BackrubSegment.cc
/// @brief implmentation of BackrubSegment class
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#include <protocols/backrub/BackrubSegment.hh>

// Core Headers
#include <core/id/AtomID.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Protocols Headers

// Utility Headers
#include <utility/keys/Key3Vector.hh>

#include <utility/vector1.hh>


// C++ Headers

using namespace core;

namespace protocols {
namespace backrub {

void
protocols::backrub::BackrubSegment::start_atoms1(
	core::pose::Pose const & pose,
	core::kinematics::tree::AtomCOP & start_atom_m1,
	core::kinematics::tree::AtomCOP & start_atom,
	core::kinematics::tree::AtomCOP & start_atom_p1
) const
{
	kinematics::AtomTree const & atom_tree(pose.atom_tree());

	start_atom = atom_tree.atom(start_atomid_).get_self_ptr();
	start_atom_m1 = start_atom->stub_atom2();
	start_atom_p1 = atom_tree.atom(start_atomid1_).get_self_ptr();
}

void
protocols::backrub::BackrubSegment::start_atoms2(
	core::pose::Pose const & pose,
	core::kinematics::tree::AtomCOP & start_atom_m2,
	core::kinematics::tree::AtomCOP & start_atom_m1,
	core::kinematics::tree::AtomCOP & start_atom,
	core::kinematics::tree::AtomCOP & start_atom_p1,
	core::kinematics::tree::AtomCOP & start_atom_p2
) const
{
	start_atoms1(pose, start_atom_m1, start_atom, start_atom_p1);

	start_atom_m2 = start_atom->stub_atom3();
	start_atom_p2 = pose.atom_tree().atom(start_atomid2_).get_self_ptr();
}

protocols::backrub::BackrubSegment::BondAngleKey
protocols::backrub::BackrubSegment::start_bond_angle_key(
	core::pose::Pose const & pose
)
{
	kinematics::tree::AtomCOP start_atom_m1, start_atom, start_atom_p1;
	start_atoms1(pose, start_atom_m1, start_atom, start_atom_p1);

	if ( start_atom_m1->id() < start_atom_p1->id() ) {
		return BondAngleKey(start_atom_m1->id(), start_atom->id(), start_atom_p1->id());
	}

	return BondAngleKey(start_atom_p1->id(), start_atom->id(), start_atom_m1->id());
}

void
protocols::backrub::BackrubSegment::end_atoms1(
	core::pose::Pose const & pose,
	core::kinematics::tree::AtomCOP & end_atom_m1,
	core::kinematics::tree::AtomCOP & end_atom,
	core::kinematics::tree::AtomCOP & end_atom_p1
) const
{
	end_atom = pose.atom_tree().atom(end_atomid_).get_self_ptr();
	end_atom_m1 = end_atom->parent();
	end_atom_p1 = end_atom->get_nonjump_atom(0);
}

void
protocols::backrub::BackrubSegment::end_atoms2(
	core::pose::Pose const & pose,
	core::kinematics::tree::AtomCOP & end_atom_m2,
	core::kinematics::tree::AtomCOP & end_atom_m1,
	core::kinematics::tree::AtomCOP & end_atom,
	core::kinematics::tree::AtomCOP & end_atom_p1,
	core::kinematics::tree::AtomCOP & end_atom_p2
) const
{
	end_atoms1(pose, end_atom_m1, end_atom, end_atom_p1);

	runtime_assert(end_atom_m1 != 0);

	end_atom_m2 = end_atom_m1->parent();
	end_atom_p2 = end_atom_p1 ? end_atom_p1->get_nonjump_atom(0) : core::kinematics::tree::AtomCOP(0) ;
}

protocols::backrub::BackrubSegment::BondAngleKey
protocols::backrub::BackrubSegment::end_bond_angle_key(
	core::pose::Pose const & pose
)
{
	kinematics::tree::AtomCOP end_atom_m1, end_atom, end_atom_p1;
	end_atoms1(pose, end_atom_m1, end_atom, end_atom_p1);

	if ( !end_atom_p1 ) {
		return BondAngleKey(end_atom_m1->id(), end_atom->id(), id::BOGUS_ATOM_ID);
	}

	if ( end_atom_m1->id() < end_atom_p1->id() ) {
		return BondAngleKey(end_atom_m1->id(), end_atom->id(), end_atom_p1->id());
	}

	return BondAngleKey(end_atom_p1->id(), end_atom->id(), end_atom_m1->id());
}

void
protocols::backrub::BackrubSegment::bond_angle_atoms(
	core::pose::Pose const & pose,
	BackrubSegment::BondAngleKey bond_angle_key,
	core::kinematics::tree::AtomCOP & atom_m1,
	core::kinematics::tree::AtomCOP & atom,
	core::kinematics::tree::AtomCOP & atom_p1
)
{
	kinematics::AtomTree const & atom_tree(pose.atom_tree());

	runtime_assert(bond_angle_key.key1().valid());
	runtime_assert(bond_angle_key.key2().valid());

	atom_m1 = atom_tree.atom(bond_angle_key.key1()).get_self_ptr();
	atom = atom_tree.atom(bond_angle_key.key2()).get_self_ptr();
	atom_p1 = bond_angle_key.key3().valid() ? atom_tree.atom(bond_angle_key.key3()).get_self_ptr() : core::kinematics::tree::AtomCOP();
}

} // moves
} // protocols
