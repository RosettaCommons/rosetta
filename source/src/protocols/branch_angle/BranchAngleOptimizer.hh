// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/branch_angle/BranchAgnleOptimizer.hh
/// @brief definition of BranchAgnleOptimizer class and methods
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_branch_angle_BranchAngleOptimizer_hh
#define INCLUDED_protocols_branch_angle_BranchAngleOptimizer_hh

#include <protocols/branch_angle/BranchAngleOptimizer.fwd.hh>

// Protocols Headers
#include <protocols/branch_angle/BranchCoef1.fwd.hh>
#include <protocols/branch_angle/BranchCoef2.fwd.hh>
#include <protocols/branch_angle/BranchParam1.fwd.hh>
#include <protocols/branch_angle/BranchParam2.fwd.hh>

#ifdef WIN32
#include <protocols/branch_angle/BranchCoef1.hh> // WIN32 INCLUDE
#include <protocols/branch_angle/BranchCoef2.hh> // WIN32 INCLUDE
#include <protocols/branch_angle/BranchParam1.hh> // WIN32 INCLUDE
#include <protocols/branch_angle/BranchParam2.hh> // WIN32 INCLUDE
#endif

// Core Headers
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers

// Standard Library Headers
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace branch_angle {

class BranchAngleOptimizer {

public:

	BranchAngleOptimizer(
		core::scoring::mm::MMBondAngleLibrary const & mm_bondangle_library
	);

	BranchAngleOptimizer();

	BranchAngleOptimizer(
		BranchAngleOptimizer const & src
	);

	core::Real
	tolerance() const
	{
		return tolerance_;
	}

	void
	tolerance(
		core::Real tolerance
	)
	{
		tolerance_ = tolerance;
	}

	virtual
	~BranchAngleOptimizer();

	core::scoring::mm::MMBondAngleResidueTypeParamSetOP
	bond_angle_residue_type_param_set();

	core::scoring::mm::MMBondAngleResidueTypeParamSetCOP
	bond_angle_residue_type_param_set() const;

	/// @brief set input MMBondAngleResidueTypeParamSet, sharing with the input
	void
	bond_angle_residue_type_param_set(
		core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set
	);

	/// @brief set input MMBondAngleResidueTypeParamSet, making a copy
	void
	bond_angle_residue_type_param_set(
		core::scoring::mm::MMBondAngleResidueTypeParamSetCOP param_set
	);

	bool initialized() const { return initialized_; }

	/// @brief optimize angles branching off the defined mainchain
	core::Size
	optimize_angles(
		core::pose::Pose & pose,
		core::id::AtomID main_atomid1,
		core::id::AtomID center_atomid,
		core::id::AtomID main_atomid2,
		bool optimize_for_minimum = false
	);

	/// @brief get overall bond angle parameters for the defined mainchian
	core::Size
	overall_params(
		core::pose::Pose const & pose,
		core::id::AtomID main_atomid1,
		core::id::AtomID center_atomid,
		core::id::AtomID main_atomid2,
		core::Real & Ktheta,
		core::Real & theta0,
		core::Real & energy0
	);

	/// @brief get single branching atom bond angle parameters
	BranchParam1
	param1(
		core::pose::Pose const & pose,
		core::id::AtomID const & main_atomid1,
		core::id::AtomID const & center_atomid,
		core::id::AtomID const & main_atomid2,
		core::id::AtomID const & branch_atomid1
	) const;

	/// @brief get double branching atom bond angle parameters
	BranchParam2
	param2(
		core::pose::Pose const & pose,
		core::id::AtomID const & main_atomid1,
		core::id::AtomID const & center_atomid,
		core::id::AtomID const & main_atomid2,
		core::id::AtomID const & branch_atomid1,
		core::id::AtomID const & branch_atomid2
	) const;

	/// @brief get number of single branching atom coefficients
	core::Size
	num_coef1() const;

	/// @brief get number of double branching atom coefficients
	core::Size
	num_coef2() const;

	/// @brief get number of undefined single branching atom coefficients
	core::Size
	num_undefined_coef1() const;

	/// @brief get number of undefined double branching atom coefficients
	core::Size
	num_undefined_coef2() const;

	/// @brief read known parameters from the database
	void
	read_database();

	/// @brief write undefined parameters to the database
	void
	write_database() const;

	/// @brief read single branching atom coefficients from an input stream
	void
	read_coef1(
		std::istream & in
	);

	/// @brief read single branching atom coefficients from a file
	bool
	read_coef1(
		std::string const & filename
	);

	/// @brief read single branching atom coefficients from default database file
	void
	read_coef1_default();

	/// @brief read single branching atom coefficients from user database file
	bool
	read_coef1_user();

	/// @brief read single branching atom coefficients
	void
	read_coef2(
		std::istream & in
	);

	/// @brief read double branching atom coefficients from a file
	bool
	read_coef2(
		std::string const & filename
	);

	/// @brief read double branching atom coefficients from default database file
	void
	read_coef2_default();

	/// @brief read double branching atom coefficients from user database file
	bool
	read_coef2_user();

	/// @brief read parameters for undefined single branching atom coefficients
	void
	read_undefined_coef1(
		std::istream & in
	);

	/// @brief read parameters for undefined single branching atom coefficients from a file
	bool
	read_undefined_coef1(
		std::string const & filename
	);

	/// @brief read single branching atom undefined coefficients from the database file
	bool
	read_undefined_coef1();

	/// @brief write out parameters for undefined single branching atom coefficients
	void
	write_undefined_coef1(
		std::ostream & out
	) const;

	/// @brief write parameters for undefined single branching atom coefficients to a file
	bool
	write_undefined_coef1(
		std::string const & filename
	) const;

	/// @brief write single branching atom undefined coefficients to the database file
	bool
	write_undefined_coef1() const;

	/// @brief read parameters for undefined double branching atom coefficients
	void
	read_undefined_coef2(
		std::istream & in
	);

	/// @brief read parameters for undefined double branching atom coefficients from a file
	bool
	read_undefined_coef2(
		std::string const & filename
	);

	/// @brief read double branching atom undefined coefficients from the database file
	bool
	read_undefined_coef2();

	/// @brief write out parameters for undefined double branching atom coefficients
	void
	write_undefined_coef2(
		std::ostream & out
	) const;

	/// @brief write parameters for undefined double branching atom coefficients to a file
	bool
	write_undefined_coef2(
		std::string const & filename
	) const;

	/// @brief write double branching atom undefined coefficients to the database file
	bool
	write_undefined_coef2() const;

private:

	// don't allow assignment
	BranchAngleOptimizer const &
	operator=(BranchAngleOptimizer const & src);

	core::scoring::mm::MMBondAngleLibrary const & mm_bondangle_library_;
	core::scoring::mm::MMBondAngleResidueTypeParamSetOP bond_angle_residue_type_param_set_;

	utility::vector1<BranchCoef1> coef1_;
	utility::vector1<BranchCoef2> coef2_;
	std::map<BranchParam1, core::Size> coef_map1_;
	std::map<BranchParam2, core::Size> coef_map2_;
	std::set<BranchParam1> undefined_coef1_;
	std::set<BranchParam2> undefined_coef2_;

	core::Real tolerance_;
	bool initialized_;
};

// Undefined, commenting out to make PyRosetta compile
/// @brief get all atoms bonded to another
//utility::vector1<core::id::AtomID> bonded_neighbor_all_res(core::pose::Pose const & pose, core::id::AtomID atomid );

/// @brief get 1 branching atom
void
branching_atomid1(
	core::pose::Pose const & pose,
	core::id::AtomID main_atomid1,
	core::id::AtomID center_atomid,
	core::id::AtomID main_atomid2,
	core::id::AtomID & branch_atomid1
);

/// @brief get 2 branching atoms ordered according their torsion offsets
void
branching_atomids2(
	core::pose::Pose const & pose,
	core::id::AtomID main_atomid1,
	core::id::AtomID center_atomid,
	core::id::AtomID main_atomid2,
	core::id::AtomID & branch_atomid1,
	core::id::AtomID & branch_atomid2
);

/// @brief get 2 siblings of an atom ordered according their torsion offsets
void
get_branching_atoms2(
	core::kinematics::tree::AtomCOP const main_atom2,
	core::kinematics::tree::AtomCOP & branch_atom1,
	core::kinematics::tree::AtomCOP & branch_atom2
);

} // branch_angle
} // protocols

#endif
