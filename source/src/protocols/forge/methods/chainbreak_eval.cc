// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/methods/chainbreak_eval.cc
/// @brief  methods for chainbreak evaluation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/methods/chainbreak_eval.hh>

// package headers
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>

// project headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/LinearChainbreakEnergy.hh>
#include <core/scoring/methods/ChainbreakEnergy.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace methods {


/// @brief evaluate linear chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
linear_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;

	Pose scratch = pose;

	// topology setup, do it here so the non-const linear_chainbreak()
	// doesn't copy the Pose a second time
	if ( !scratch.fold_tree().is_cutpoint( pos ) || scratch.fold_tree().num_cutpoint() > 1 ) {
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );
	}

	return linear_chainbreak( scratch, pos ); // call non-const
}


/// @brief evaluate linear chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.  If cutpoint variants are present at chainbreak, will use
///  existing variants and not modify them.  If cutpoint variants are not
///  found will add them and then remove them once calculation is finished.
core::Real
linear_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;
	using core::scoring::EnergyMap;
	using core::scoring::ScoreFunction;
	using core::scoring::ScoreType;
	using core::scoring::methods::LinearChainbreakEnergy;

	if ( pose.fold_tree().num_cutpoint() == 0 ) {
		return 0;
	}
	assert( pos > 0 );
	assert( pos < pose.n_residue() );

	EnergyMap emap;
	ScoreFunction fx; // dummy, needed for function call
	LinearChainbreakEnergy energy;

	// evaluate the break
	if ( pose.fold_tree().is_cutpoint( pos ) && pose.fold_tree().num_cutpoint() == 1 ) {

		bool cutpoints_added = add_cutpoint_variants( pose, pos );
		energy.finalize_total_energy( pose, fx, emap );
		if ( cutpoints_added ) {
			remove_cutpoint_variants( pose, pos );
		}

	} else { // no cutpoint, copy the Pose and introduce one

		Pose scratch = pose;
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );

		add_cutpoint_variants( scratch, pos );
		energy.finalize_total_energy( scratch, fx, emap );

	}

	// return the energy
	return emap[ core::scoring::linear_chainbreak ];
}


/// @brief evaluate overlap chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
overlap_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;

	Pose scratch = pose;

	// topology setup, do it here so the non-const linear_chainbreak()
	// doesn't copy the Pose a second time
	if ( !scratch.fold_tree().is_cutpoint( pos ) || scratch.fold_tree().num_cutpoint() > 1 ) {
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );
	}

	return overlap_chainbreak( scratch, pos ); // call non-const
}


/// @brief evaluate overlap chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.  If cutpoint variants are present at chainbreak, will use
///  existing variants and not modify them.  If cutpoint variants are not
///  found will add them and then remove them once calculation is finished.
core::Real
overlap_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;
	using core::scoring::EnergyMap;
	using core::scoring::ScoreFunction;
	using core::scoring::ScoreType;
	using core::scoring::methods::LinearChainbreakEnergy;

	assert( pos > 0 );
	assert( pos < pose.n_residue() );

	EnergyMap emap;
	ScoreFunction fx; // dummy, needed for function call
	LinearChainbreakEnergy energy; // also calculates the overlap chainbreak

	// evaluate the break
	if ( pose.fold_tree().is_cutpoint( pos ) && pose.fold_tree().num_cutpoint() == 1 ) {

		bool cutpoints_added = add_cutpoint_variants( pose, pos );
		energy.finalize_total_energy( pose, fx, emap );
		if ( cutpoints_added ) {
			remove_cutpoint_variants( pose, pos );
		}

	} else { // no cutpoint, copy the Pose and introduce one

		Pose scratch = pose;
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );

		add_cutpoint_variants( scratch, pos );
		energy.finalize_total_energy( scratch, fx, emap );

	}

	// return the energy
	return emap[ core::scoring::overlap_chainbreak ];
}


/// @brief evaluate quadratic chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
quadratic_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;

	Pose scratch = pose;

	// topology setup, do it here so the non-const quadratic_chainbreak()
	// doesn't copy the Pose a second time
	if ( !scratch.fold_tree().is_cutpoint( pos ) || scratch.fold_tree().num_cutpoint() > 1 ) {
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );
	}

	return quadratic_chainbreak( scratch, pos ); // call non-const
}


/// @brief evaluate quadratic chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.  If cutpoint variants are present at chainbreak, will use
///  existing variants and not modify them.  If cutpoint variants are not
///  found will add them and then remove them once calculation is finished.
core::Real
quadratic_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
)
{
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;
	using core::scoring::EnergyMap;
	using core::scoring::ScoreFunction;
	using core::scoring::ScoreType;
	using core::scoring::methods::ChainbreakEnergy;

	assert( pos > 0 );
	assert( pos < pose.n_residue() );

	EnergyMap emap;
	ScoreFunction fx; // dummy, needed for function call
	ChainbreakEnergy energy;

	// evaluate the break
	if ( pose.fold_tree().is_cutpoint( pos ) && pose.fold_tree().num_cutpoint() == 1 ) {

		bool cutpoints_added = add_cutpoint_variants( pose, pos );
		energy.finalize_total_energy( pose, fx, emap );
		if ( cutpoints_added ) {
			remove_cutpoint_variants( pose, pos );
		}

	} else { // no cutpoint, copy the Pose and introduce one

		Pose scratch = pose;
		FoldTree ft = fold_tree_from_pose( scratch, scratch.fold_tree().root(), MoveMap() );
		ft.new_jump( pos, pos + 1, pos );
		scratch.fold_tree( ft );

		add_cutpoint_variants( scratch, pos );
		energy.finalize_total_energy( scratch, fx, emap );

	}

	// return the energy
	return emap[ core::scoring::chainbreak ];
}


} // methods
} // forge
} // protocols
