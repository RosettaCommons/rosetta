// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/branch_angle/BranchParam1.hh
/// @brief definition/implementation of BranchParam1 class and methods
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_branch_angle_BranchParam1_hh
#define INCLUDED_protocols_branch_angle_BranchParam1_hh

#include <protocols/branch_angle/BranchParam1.fwd.hh>

// Core Headers
#include <core/types.hh>

namespace protocols {
namespace branch_angle {

/// @brief
/// a class to store bond angle energy parameters around a single atom atom with three bonded neighbors
class BranchParam1 {

public:

	BranchParam1(
		core::Real m1_m2_Ktheta,
		core::Real m1_m2_theta0,
		core::Real m1_b1_Ktheta,
		core::Real m1_b1_theta0,
		core::Real m2_b1_Ktheta,
		core::Real m2_b1_theta0,
		core::Real tolerance = 0
	):
		m1_m2_Ktheta_(m1_m2_Ktheta),
		m1_m2_theta0_(m1_m2_theta0),
		m1_b1_Ktheta_(m1_b1_Ktheta),
		m1_b1_theta0_(m1_b1_theta0),
		m2_b1_Ktheta_(m2_b1_Ktheta),
		m2_b1_theta0_(m2_b1_theta0),
		tolerance_(tolerance)
	{}

	/// @brief get Ktheta for mainchain atom 1 - mainchain atom 1 angle
	core::Real
	m1_m2_Ktheta() const
	{
		return m1_m2_Ktheta_;
	}

	/// @brief get theta0 for mainchain atom 1 - mainchain atom 1 angle
	core::Real
	m1_m2_theta0() const
	{
		return m1_m2_theta0_;
	}

	/// @brief get Ktheta for mainchain atom 1 - branching atom 1 angle
	core::Real
	m1_b1_Ktheta() const
	{
		return m1_b1_Ktheta_;
	}

	/// @brief get theta0 for mainchain atom 1 - branching atom 1 angle
	core::Real
	m1_b1_theta0() const
	{
		return m1_b1_theta0_;
	}

	/// @brief get Ktheta for mainchain atom 2 - branching atom 1 angle
	core::Real
	m2_b1_Ktheta() const
	{
		return m2_b1_Ktheta_;
	}

	/// @brief get theta0 for mainchain atom 2 - branching atom 1 angle
	core::Real
	m2_b1_theta0() const
	{
		return m2_b1_theta0_;
	}

	/// @brief a is LOWER than b by a given tolerance
	friend
	inline
	bool
	operator <(
		BranchParam1 const & a,
		BranchParam1 const & b
	)
	{
		core::Real const tolerance(a.tolerance_ < b.tolerance_ ? a.tolerance_ : b.tolerance_);

		core::Real const m1_m2_Ktheta_diff(a.m1_m2_Ktheta_ - b.m1_m2_Ktheta_);
		if ( m1_m2_Ktheta_diff < -tolerance ) return true;
		if ( m1_m2_Ktheta_diff > tolerance ) return false;
		// within tolerance: proceed to checking the next number

		core::Real const m1_m2_theta0_diff(a.m1_m2_theta0_ - b.m1_m2_theta0_);
		if ( m1_m2_theta0_diff < -tolerance ) return true;
		if ( m1_m2_theta0_diff > tolerance ) return false;
		// within tolerance: proceed to checking the next number

		core::Real const m1_b1_Ktheta_diff(a.m1_b1_Ktheta_ - b.m1_b1_Ktheta_);
		if ( m1_b1_Ktheta_diff < -tolerance ) return true;
		if ( m1_b1_Ktheta_diff > tolerance ) return false;
		// within tolerance: proceed to checking the next number

		core::Real const m1_b1_theta0_diff(a.m1_b1_theta0_ - b.m1_b1_theta0_);
		if ( m1_b1_theta0_diff < -tolerance ) return true;
		if ( m1_b1_theta0_diff > tolerance ) return false;
		// within tolerance: proceed to checking the next number

		core::Real const m2_b1_Ktheta_diff(a.m2_b1_Ktheta_ - b.m2_b1_Ktheta_);
		if ( m2_b1_Ktheta_diff < -tolerance ) return true;
		if ( m2_b1_Ktheta_diff > tolerance ) return false;
		// within tolerance: proceed to checking the next number

		core::Real const m2_b1_theta0_diff(a.m2_b1_theta0_ - b.m2_b1_theta0_);
		if ( m2_b1_theta0_diff < -tolerance ) return true;
		// either greater than or within tolarance
		return false;
	}

protected:

	core::Real m1_m2_Ktheta_;
	core::Real m1_m2_theta0_;
	core::Real m1_b1_Ktheta_;
	core::Real m1_b1_theta0_;
	core::Real m2_b1_Ktheta_;
	core::Real m2_b1_theta0_;

	core::Real tolerance_;
};

} // branch_angle
} // protocols

#endif
