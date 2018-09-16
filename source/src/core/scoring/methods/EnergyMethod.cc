// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/method/EnergyMethod.cc
/// @brief  Base class for energy classes
/// @author Phil Bradley


// Unit headers
#include <core/scoring/methods/EnergyMethod.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

EnergyMethod::EnergyMethod( EnergyMethodCreatorOP creator ) :
	score_types_( creator->score_types_for_method() )
{
}

EnergyMethod::EnergyMethod( EnergyMethod const & src ) :
	parent(),
	score_types_( src.score_types_ )
{
}

EnergyMethod::~EnergyMethod() = default;

bool
EnergyMethod::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	/// APL -- Grandfathering mechanism.  As energy methods are updated to use
	/// the new derivative evaluation machinery, they can set this function to
	/// return "false", and therefore be used in rtmin.
	return true;
}

bool
EnergyMethod::defines_high_order_terms( pose::Pose const & ) const
{
	return false;
}

/// @details default implementation noop
void
EnergyMethod::setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const {}

/// @brief if an energy method needs to cache data in the Energies object,
/// before packing begins and requires access to the RotamerSets object, then
/// it does so during this function. The default behavior is to do nothing.
/// @details The exact order of events when setting up for packing are as follows:
///          1. setup_for_packing() is called for all energy methods
///          2. rotamers are built
///          3. setup_for_packing_with_rotsets() is called for all energy methods
///          4. prepare_rotamers_for_packing() is called for all energy methods
///          5. The energy methods are asked to score all rotamers and rotamer pairs
///          6. Annealing
/// @remarks The pose is specifically non-const here so that energy methods can store data in it
/// @note: Used in ApproximateBuriedUnsatPenalty to pre-compute compatible rotamers
void
EnergyMethod::setup_for_packing_with_rotsets(
	pose::Pose &,
	pack_basic::RotamerSetsBaseOP const &,
	ScoreFunction const &
) const {}

/// @details default implementation noop
void
EnergyMethod::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase & ) const
{}


/// @details default implementation noop
void
EnergyMethod::update_residue_for_packing(
	pose::Pose &,
	Size ) const
{}

void
EnergyMethod::setup_for_scoring( pose::Pose &, ScoreFunction const & ) const {}

bool
EnergyMethod::requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose::Pose const & ) const
{
	return false;
}

void
EnergyMethod::setup_for_scoring_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	basic::datacache::BasicDataCache &
) const
{}


/// @details default implementation noop
void
EnergyMethod::setup_for_minimizing(
	pose::Pose & ,
	ScoreFunction const & ,
	kinematics::MinimizerMapBase const &
) const {}

/// @brief Called after minimization, allowing a derived class to do some
/// teardown steps.
/// @details Base class function does nothing.  Derived classes may override.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
EnergyMethod::finalize_after_minimizing(
	pose::Pose &// pose
) const {}

/// @details default implementation noop
void
EnergyMethod::setup_for_derivatives( pose::Pose &, ScoreFunction const &  ) const {}

/// @details default implementation noop
void
EnergyMethod::finalize_after_derivatives( pose::Pose &, ScoreFunction const &  ) const
{}

/// @details default implementation does not alter either the F1 or F2 vectors.
void
EnergyMethod::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,// F1,
	Vector & // F2
) const {}

/// called at the end of energy evaluation
void
EnergyMethod::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const {}

/// @brief Override the entirety of the score types list if they
/// were initialized incorrectly in a parent's constructor.
void
EnergyMethod::set_score_types( EnergyMethodCreatorOP creator ) {
	score_types_ = creator->score_types_for_method();
}

/// @brief show additional information of the energy method
void EnergyMethod::show_additional_info(std::ostream &, pose::Pose &, bool) const {}

} // methods
} // scoring
} // core
