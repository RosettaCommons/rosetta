// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <core/scoring/methods/OneBodyEnergy.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.OneBodyEnergy" );

OneBodyEnergy::OneBodyEnergy(
	EnergyMethodCreatorOP creator
) :
	parent( creator )
{}

OneBodyEnergy::~OneBodyEnergy() {}

bool
OneBodyEnergy::defines_score_for_residue(
	conformation::Residue const &
) const
{
	return true;
}

bool
OneBodyEnergy::use_extended_residue_energy_interface() const
{
	return false;
}

void
OneBodyEnergy::residue_energy_ext(
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap &
) const
{}


void
OneBodyEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const &,
	pose::Pose const & ,
	ScoreFunction const & ,
	kinematics::MinimizerMapBase const & ,
	ResSingleMinimizationData &
) const {}

bool
OneBodyEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose::Pose const & ) const
{
	return false;
}

void
OneBodyEnergy::setup_for_scoring_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData &
) const
{
	// noop -- this should be an error
}

bool
OneBodyEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & ) const
{
	return false;
}

void
OneBodyEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData &,
	basic::datacache::BasicDataCache &
) const
{
	// noop -- this should be an error since this function could only be called if the
	// derived class's requires_a_setup_for_derivatives_for_residue_opportunity message
	// returns true, but the derived class has not implemented this function or has
	// improperly overriden the base class version.
	TR << "WARNING: Unimplemented or improperly overridden OneBodyEnergy::setup_for_derivatives_for_residue() for class computing";
	for ( auto st : score_types() ) {
		TR << " " << st;
	}
	TR << std::endl;
}


void
OneBodyEnergy::eval_residue_derivatives(
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const &,
	utility::vector1< DerivVectorPair > &
) const
{}

bool
OneBodyEnergy::defines_dof_derivatives( pose::Pose const & ) const
{
	return false;
}

Real
OneBodyEnergy::eval_residue_dof_derivative(
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}


} // methods
} // scoring
} // core

