// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/ConstantConstraint.hh

#ifndef INCLUDED_core_scoring_constraints_ConstantConstraint_hh
#define INCLUDED_core_scoring_constraints_ConstantConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

/// @brief A Constant Constraint.
class ConstantConstraint : public Constraint {
public:
	using Constraint::score;

	/// @brief compute score
	Real
	score() const;

	virtual std::string type() const;

	/// @brief compute score
	void
	score( func::XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const;

	//Real score( pose::Pose const & ) const;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		AtomID const & ,
		func::XYZ_Func const &,
		Vector & F1,
		Vector & F2,
		EnergyMap const &
	) const;

	/// @brief Constructor
	ConstantConstraint(
		func::FuncOP func_in,
		ScoreType scotype = constant_constraint
	);

	// destructor
	~ConstantConstraint();

	/// @brief number of atoms --- zero
	Size
	natoms() const;

	AtomID const &
	atom( Size const n ) const;

	/// @brief output violation of constraint (none!)
	Size show_violations( std::ostream &, pose::Pose const &, Size, core::Real ) const;

	void show( std::ostream& out ) const;

private:
	func::FuncOP func_;
};

} // constraints
} // scoring
} // core

#endif
