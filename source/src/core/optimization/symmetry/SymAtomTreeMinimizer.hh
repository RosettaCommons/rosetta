// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/AtomTreeMinimizer.hh
/// @brief  High-level atom tree minimizer class for symmetrical minimization
/// @author Ingemar Andre

#ifndef INCLUDED_core_optimization_symmetry_SymAtomTreeMinimizer_hh
#define INCLUDED_core_optimization_symmetry_SymAtomTreeMinimizer_hh


// Package headers
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
// Symmetry
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {
namespace symmetry {

/// @brief High-level atom tree minimizer class
class SymAtomTreeMinimizer : public core::optimization::AtomTreeMinimizer
{

public:

	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfo SymmetryInfoOP;

public:

	// c-tor
	SymAtomTreeMinimizer(){};

	/// @brief Override the base class implementation.  Non-const.
	virtual
	Real
	run(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map,
		scoring::ScoreFunction const & scorefxn,
		MinimizerOptions const & options
	);

	static
	void
	make_assymetric_movemap(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map_sym,
		kinematics::MoveMap & move_map_asym
	);

	static
	void
	make_semisymmetric_movemap(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map_sym,
		kinematics::MoveMap & move_map_semisym
	);

}; // SymAtomTreeMinimizer

} // symmetry
} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_symmetry_AtomTreeMinimizer_HH
