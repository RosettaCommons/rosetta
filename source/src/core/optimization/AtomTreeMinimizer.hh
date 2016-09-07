// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/AtomTreeMinimizer.hh
/// @brief  High-level atom tree minimizer class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_AtomTreeMinimizer_hh
#define INCLUDED_core_optimization_AtomTreeMinimizer_hh

// Unit headers
#include <core/optimization/AtomTreeMinimizer.fwd.hh>

// Package headers
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {


/// @brief High-level atom tree minimizer class
class AtomTreeMinimizer : public utility::pointer::ReferenceCount
{

public:

	// c-tor
	AtomTreeMinimizer();

	~AtomTreeMinimizer() override;

	/// @brief run minimization and return the final score at minimization's conclusion.
	/// Virtual allowing derived classes to mascarade as AtomTreeMinimizers.
	/// Non-const so that it can modify its deriv_check_result_ object.
	virtual
	Real
	run(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map,
		scoring::ScoreFunction const & scorefxn,
		MinimizerOptions const & options
	);

	/// @brief After minimization has concluded, the user may access the deriv-check result,
	/// assuming that they have run the AtomTreeMinimizer with deriv_check = true;
	NumericalDerivCheckResultOP
	deriv_check_result() const;

	/// @brief Do consistency checks for minimizer setup.
	void
	check_setup(pose::Pose const & pose,
		kinematics::MoveMap const & move_map,
		scoring::ScoreFunction const & scorefxn,
		MinimizerOptions const & options) const;

private:

	NumericalDerivCheckResultOP deriv_check_result_;

}; // AtomTreeMinimizer


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_AtomTreeMinimizer_HH
