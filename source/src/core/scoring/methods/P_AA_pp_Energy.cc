// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/P_AA_pp_Energy.cc
/// @brief  Probability of observing an amino acid, given its phi/psi energy method declaration
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/P_AA_pp_Energy.hh>
#include <core/scoring/methods/P_AA_pp_EnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/conversions.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
P_AA_pp_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new P_AA_pp_Energy );
}

ScoreTypes
P_AA_pp_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( p_aa_pp );
	return sts;
}


/// ctor
P_AA_pp_Energy::P_AA_pp_Energy() :
	parent( methods::EnergyMethodCreatorOP( new P_AA_pp_EnergyCreator ) ),
	p_aa_( ScoringManager::get_instance()->get_P_AA() )
{}

/// clone
EnergyMethodOP
P_AA_pp_Energy::clone() const
{
	return EnergyMethodOP( new P_AA_pp_Energy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
P_AA_pp_Energy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	emap[ p_aa_pp ] += p_aa_.P_AA_pp_energy( rsd );
}


bool
P_AA_pp_Energy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}


Real
P_AA_pp_Energy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,// min_data,
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const &,// pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	if ( ! tor_id.valid() ) return 0.0;
	return numeric::conversions::degrees( weights[ p_aa_pp ] * p_aa_.get_Paa_pp_deriv( rsd, tor_id ));
}


Real
P_AA_pp_Energy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &, //sfxn,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	if ( ! tor_id.valid() ) return 0.0;
	return numeric::conversions::degrees( weights[ p_aa_pp ] * p_aa_.get_Paa_pp_deriv( pose.residue( tor_id.rsd() ), tor_id ));
}

/// @brief P_AA_pp_Energy is context independent; indicates that no
/// context graphs are required
void
P_AA_pp_Energy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
P_AA_pp_Energy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

