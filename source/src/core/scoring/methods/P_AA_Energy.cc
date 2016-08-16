// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/P_AA_Energy.cc
/// @brief  Probability of observing an amino acid (NOT conditional on phi/psi), energy method implementation
/// @author Ron Jacak


// Unit headers
#include <core/scoring/methods/P_AA_Energy.hh>
#include <core/scoring/methods/P_AA_EnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <utility/vector1.hh>


// Project headers


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the P_AA_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
P_AA_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new P_AA_Energy );
}

ScoreTypes
P_AA_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( p_aa );
	return sts;
}


/// @remarks
/// get_P_AA calls a method in ScoringManager which create a new object of type P_AA.  The constructor for that created object
/// reads in the three database files: p_aa, p_aa_pp, and p_aa_n.  That object is returned and then stored as a private member
/// variable here.
P_AA_Energy::P_AA_Energy() :
	parent( methods::EnergyMethodCreatorOP( new P_AA_EnergyCreator ) ),
	p_aa_( ScoringManager::get_instance()->get_P_AA() )
{}


EnergyMethodOP
P_AA_Energy::clone() const {
	return EnergyMethodOP( new P_AA_Energy );
}


/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

void
P_AA_Energy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	emap[ p_aa ] += p_aa_.P_AA_energy( rsd );
}


/// @remarks no DOF to vary for P_AA, so just return 0.0 like the reference energy does
Real
P_AA_Energy::eval_dof_derivative(
	id::DOF_ID const & ,// dof_id,
	id::TorsionID const & , //tor_id
	pose::Pose const & , // pose
	ScoreFunction const &, //sfxn,
	EnergyMap const & // weights
) const
{
	return 0.0;
}

/// @brief P_AA_Energy is context independent; indicates that no context graphs are required
void
P_AA_Energy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}
core::Size
P_AA_Energy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

