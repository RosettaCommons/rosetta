// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/OmegaTetherEnergy.hh
/// @brief  OmegaTether energy method class declaration
/// @details  This score term constrains the inter-residue torsion (omega) to be 0 or 180 degrees.
/// It works for alpha-amino acids, beta-amino acids, and oligoureas.  In the case of oligoureas,
/// it constrains both omega and mu (the preceding torsion) to be 180.
/// @author Phil Bradley


#ifndef INCLUDED_core_energy_methods_OmegaTetherEnergy_hh
#define INCLUDED_core_energy_methods_OmegaTetherEnergy_hh

// Unit headers
#include <core/energy_methods/OmegaTetherEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class OmegaTetherEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;
public:

	/// ctor
	OmegaTetherEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	/// @brief Use the dof_derivative interface for this energy method when
	/// calculating derivatives?  It is possible to define both dof_derivatives and
	/// atom-derivatives; they are not mutually exclusive.
	bool
	defines_dof_derivatives( pose::Pose const & p ) const override;

	/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
	/// and the input residue is not required to be a member of the Pose.
	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const override;

	/// @brief OmegaTether Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	// data
private:
	core::scoring::OmegaTether const & potential_;
	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
