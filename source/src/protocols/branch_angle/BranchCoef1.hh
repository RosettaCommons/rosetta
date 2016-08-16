// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/branch_angle/BranchCoef1.hh
/// @brief definition/implementation of BranchCoef1 class and methods
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_branch_angle_BranchCoef1_hh
#define INCLUDED_protocols_branch_angle_BranchCoef1_hh

#include <protocols/branch_angle/BranchCoef1.fwd.hh>

// Core Headers
#include <core/types.hh>

namespace protocols {
namespace branch_angle {

/// @brief
/// a class to store coefficients for branching angle optimization around a
/// single atom atom with three bonded neighbors
class BranchCoef1 {

public:

	BranchCoef1(
		core::Real overall_Ktheta,
		core::Real overall_theta0,
		core::Real overall_energy0,
		core::Real b1_torsion_offset_A,
		core::Real b1_torsion_offset_B,
		core::Real b1_torsion_offset_C,
		core::Real b1_bond_angle_A,
		core::Real b1_bond_angle_B,
		core::Real b1_bond_angle_C
	):
		overall_Ktheta_(overall_Ktheta),
		overall_theta0_(overall_theta0),
		overall_energy0_(overall_energy0),
		b1_torsion_offset_A_(b1_torsion_offset_A),
		b1_torsion_offset_B_(b1_torsion_offset_B),
		b1_torsion_offset_C_(b1_torsion_offset_C),
		b1_bond_angle_A_(b1_bond_angle_A),
		b1_bond_angle_B_(b1_bond_angle_B),
		b1_bond_angle_C_(b1_bond_angle_C)
	{}

	/// @brief get overall Ktheta parameter
	core::Real
	overall_Ktheta() const
	{
		return overall_Ktheta_;
	}

	/// @brief get overall theta0 parameter
	core::Real
	overall_theta0() const
	{
		return overall_theta0_;
	}

	/// @brief get overall energy0 parameter
	core::Real
	overall_energy0() const
	{
		return overall_energy0_;
	}

	/// @brief get branching atom 1 torsion offset A coefficient (angle^0)
	core::Real
	b1_torsion_offset_A() const
	{
		return b1_torsion_offset_A_;
	}

	/// @brief get branching atom 1 torsion offset B coefficient (angle^1)
	core::Real
	b1_torsion_offset_B() const
	{
		return b1_torsion_offset_B_;
	}

	/// @brief get branching atom 1 torsion offset C coefficient (angle^2)
	core::Real
	b1_torsion_offset_C() const
	{
		return b1_torsion_offset_C_;
	}

	/// @brief get branching atom 1 bond angle A coefficient (angle^0)
	core::Real
	b1_bond_angle_A() const
	{
		return b1_bond_angle_A_;
	}

	/// @brief get branching atom 1 bond angle B coefficient (angle^1)
	core::Real
	b1_bond_angle_B() const
	{
		return b1_bond_angle_B_;
	}

	/// @brief get branching atom 1 bond angle C coefficient (angle^2)
	core::Real
	b1_bond_angle_C() const
	{
		return b1_bond_angle_C_;
	}

	/// @brief calculate single branching angles for a main chain bond angle
	void
	evaluate(
		core::Real const m2_bond_angle,
		core::Real & b1_torsion_offset,
		core::Real & b1_bond_angle
	) const
	{
		b1_torsion_offset = b1_torsion_offset_A_ +
			m2_bond_angle * (b1_torsion_offset_B_ + m2_bond_angle * b1_torsion_offset_C_);
		b1_bond_angle = b1_bond_angle_A_ + m2_bond_angle * (b1_bond_angle_B_ + m2_bond_angle * b1_bond_angle_C_);
	}

private:

	core::Real overall_Ktheta_;
	core::Real overall_theta0_;
	core::Real overall_energy0_;
	core::Real b1_torsion_offset_A_;
	core::Real b1_torsion_offset_B_;
	core::Real b1_torsion_offset_C_;
	core::Real b1_bond_angle_A_;
	core::Real b1_bond_angle_B_;
	core::Real b1_bond_angle_C_;
};

} // branch_angle
} // protocols

#endif
