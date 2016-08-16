// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backrub/BackrubSegment.hh
/// @brief definition/implmentation of BackrubSegment class
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_backrub_BackrubSegment_hh
#define INCLUDED_protocols_backrub_BackrubSegment_hh

// Core Headers
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

// Protocols Headers

// Utility Headers
#include <utility/keys/Key3Vector.fwd.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace protocols {
namespace backrub {

/// @brief a class for holind information about individual backrub segments
class BackrubSegment {

public:

	typedef utility::keys::Key3Vector<core::id::AtomID> BondAngleKey;

	BackrubSegment(
		core::id::AtomID start_atomid,
		core::id::AtomID start_atomid1,
		core::id::AtomID start_atomid2,
		core::id::AtomID end_atomid,
		core::Size size,
		core::Real max_angle_disp
	):
		start_atomid_(start_atomid),
		start_atomid1_(start_atomid1),
		start_atomid2_(start_atomid2),
		end_atomid_(end_atomid),
		size_(size),
		max_angle_disp_(max_angle_disp),
		angle_disp_(0)
	{}

	/// @brief get AtomID of starting atom
	core::id::AtomID
	start_atomid() const
	{
		return start_atomid_;
	}

	/// @brief get AtomID of first atom along the path from start to end
	core::id::AtomID
	start_atomid1() const
	{
		return start_atomid1_;
	}

	/// @brief get AtomID of second atom along the path from start to end
	core::id::AtomID
	start_atomid2() const
	{
		return start_atomid2_;
	}

	/// @brief get AtomID of starting atom
	core::id::AtomID
	end_atomid() const
	{
		return end_atomid_;
	}

	/// @brief get the segment size
	core::Size
	size() const
	{
		return size_;
	}

	/// @brief get mainchain atom tree atoms 1 bond away from the start pivot
	void
	start_atoms1(
		core::pose::Pose const & pose,
		core::kinematics::tree::AtomCOP & start_atom_m1,
		core::kinematics::tree::AtomCOP & start_atom,
		core::kinematics::tree::AtomCOP & start_atom_p1
	) const;

	/// @brief get mainchain atom tree atoms 2 bonds away from the start pivot
	void
	start_atoms2(
		core::pose::Pose const & pose,
		core::kinematics::tree::AtomCOP & start_atom_m2,
		core::kinematics::tree::AtomCOP & start_atom_m1,
		core::kinematics::tree::AtomCOP & start_atom,
		core::kinematics::tree::AtomCOP & start_atom_p1,
		core::kinematics::tree::AtomCOP & start_atom_p2
	) const;

	/// @brief get a key representing the starting mainchain bond angle atoms
	BondAngleKey
	start_bond_angle_key(
		core::pose::Pose const & pose
	);

	/// @brief get mainchain atom tree atoms 1 bond away from the end pivot
	void
	end_atoms1(
		core::pose::Pose const & pose,
		core::kinematics::tree::AtomCOP & end_atom_m1,
		core::kinematics::tree::AtomCOP & end_atom,
		core::kinematics::tree::AtomCOP & end_atom_p1
	) const;

	/// @brief get mainchain atom tree atoms 2 bonds away from the end pivot
	void
	end_atoms2(
		core::pose::Pose const & pose,
		core::kinematics::tree::AtomCOP & end_atom_m2,
		core::kinematics::tree::AtomCOP & end_atom_m1,
		core::kinematics::tree::AtomCOP & end_atom,
		core::kinematics::tree::AtomCOP & end_atom_p1,
		core::kinematics::tree::AtomCOP & end_atom_p2
	) const;

	/// @brief get a key representing the ending mainchain bond angle atoms
	BondAngleKey
	end_bond_angle_key(
		core::pose::Pose const & pose
	);

	/// @brief get the current bond angle atoms referred to by a key
	static
	void
	bond_angle_atoms(
		core::pose::Pose const & pose,
		BackrubSegment::BondAngleKey bond_angle_key,
		core::kinematics::tree::AtomCOP & atom_m1,
		core::kinematics::tree::AtomCOP & atom,
		core::kinematics::tree::AtomCOP & atom_p1
	);

	/// @brief get maximum angular displacement
	core::Real
	max_angle_disp() const
	{
		return max_angle_disp_;
	}

	/// @brief get overall angular displacement
	core::Real
	angle_disp() const
	{
		return angle_disp_;
	}

private:

	core::id::AtomID start_atomid_;
	core::id::AtomID start_atomid1_;
	core::id::AtomID start_atomid2_;
	core::id::AtomID end_atomid_;
	core::Size size_;
	core::Real max_angle_disp_;
	core::Real angle_disp_;
	// to be implemented
	//utility::histogram<float> > attempted_moves_;
	//utility::histogram<float> > accepted_moves_;
};

} // moves
} // protocols

#endif
