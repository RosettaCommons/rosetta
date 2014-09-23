// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/OmegaTetherEnergy.cc
/// @brief  OmegaTether energy method class implementation
/// @author Phil Bradley
/// @author Mike Tyka (mtyka@u.washington.edu)

// Unit Headers
#include <core/scoring/methods/OmegaTetherEnergy.hh>
#include <core/scoring/methods/OmegaTetherEnergyCreator.hh>

// Package Headers
#include <core/scoring/OmegaTether.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>


// Utility headers
#include <numeric/conversions.hh>

#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>


// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the OmegaTetherEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
OmegaTetherEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new OmegaTetherEnergy );
}

ScoreTypes
OmegaTetherEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( omega );
	return sts;
}



/// ctor
OmegaTetherEnergy::OmegaTetherEnergy() :
	parent( methods::EnergyMethodCreatorOP( new OmegaTetherEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_OmegaTether() )
{}

/// clone
EnergyMethodOP
OmegaTetherEnergy::clone() const
{
	return EnergyMethodOP( new OmegaTetherEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
OmegaTetherEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;

	if ( rsd.is_protein() ) {
		Real omega_score, dscore_domega, dscore_dphi, dscore_dpsi;
		potential_.eval_omega_score_residue( rsd, omega_score, dscore_domega, dscore_dphi, dscore_dpsi );
		emap[ omega ] += omega_score;
	}
}


/// @brief Use the dof_derivative interface for this energy method when
/// calculating derivatives?  It is possible to define both dof_derivatives and
/// atom-derivatives; they are not mutually exclusive.
bool
OmegaTetherEnergy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}

/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
/// and the input residue is not required to be a member of the Pose.
Real
OmegaTetherEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	id::DOF_ID const &,
	id::TorsionID const & tor_id,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return 0.0;
	}

	Real deriv(0.0);
	if ( tor_id.valid() && tor_id.type() == id::BB && tor_id.torsion() <= ((rsd.has("CM") && rsd.mainchain_torsions().size()==4) ? 4 : 3)  && rsd.is_protein() ) {
		Real omega_score, dscore_domega, dscore_dphi, dscore_dpsi;
		potential_.eval_omega_score_residue( rsd, omega_score, dscore_domega, dscore_dphi, dscore_dpsi );
		if (tor_id.torsion() == 1) deriv = dscore_dphi;
		if (tor_id.torsion() == 2) deriv = dscore_dpsi;
		if ((rsd.has("CM") && rsd.mainchain_torsions().size()==4 && tor_id.torsion() == 4) || tor_id.torsion() == 3) deriv = dscore_domega; //If this is tor 3, or if this is a beta-amino acid residue and this is tor 4, store the derivative WRT omega.
	}
	return numeric::conversions::degrees( weights[ omega ] * deriv );
}

///
Real
OmegaTetherEnergy::old_eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ){
			return 0.0;
	}

	Real deriv(0.0);
	if ( tor_id.valid() && tor_id.type() == id::BB ) {
		conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
		if ( rsd.is_protein() &&
				 tor_id.torsion() <= 3 ) {
			Real omega_score, dscore_domega, dscore_dphi, dscore_dpsi;
			potential_.eval_omega_score_residue( rsd, omega_score, dscore_domega, dscore_dphi, dscore_dpsi );
			if (tor_id.torsion() == 1) deriv = dscore_dphi;
			if (tor_id.torsion() == 2) deriv = dscore_dpsi;
			if (tor_id.torsion() == 3) deriv = dscore_domega;
		}
	}
	// note that the atomtree Oomega dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ omega ] * deriv );
}

/// @brief OmegaTether Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
OmegaTetherEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
OmegaTetherEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

