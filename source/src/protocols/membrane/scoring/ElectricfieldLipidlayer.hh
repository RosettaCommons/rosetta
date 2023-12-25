// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/ElectricfieldLipidlayer.hh
/// @brief Implicit Lipid Membrane Model electrostatic energy due to the field created by lipid layers(one-body)
/// @author  Rituparna Samanta (rsamant2@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_hh
#define INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_hh

// Unit headers
#include <protocols/membrane/scoring/ElectricfieldLipidlayer.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Package headers
#include <protocols/membrane/scoring/MEnvElectroAtomParams.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief Fullatom Membrane Environment Energy
class ElectricfieldLipidlayer : public core::scoring::methods::ContextDependentOneBodyEnergy {

public:

	typedef core::scoring::methods::ContextDependentOneBodyEnergy parent;

	/// @brief Construct Energy Method from Etable
	ElectricfieldLipidlayer();

	/// @brief Clone Energy Method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief Compute Per-Residue Energies
	void
	residue_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;


	/// @brief Setup for Computing Derivatives
	void
	setup_for_derivatives(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const & scfxn
	) const override;

	/// @brief Evaluate Per-Atom Derivatives
	void
	eval_atom_derivative(
		core::id::AtomID const & id,
		core::pose::Pose const & pose,
		core::kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & emap,
		core::Vector & F1,
		core::Vector & F2
	) const override;

	/// @brief Fa_MbenvEnergy is context independent
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	/// @brief Setup Method for initial scoring
	void
	setup_for_scoring(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &
	) const override;

	/// @brief Evaluate Per-Atom Energy term
	core::Real
	eval_felec_lipidlayer(
		MEnvElectroAtomParams const & p
	) const;
	/// @brief Evaluate Per-Atom Env term
	MEnvElectroAtomParamsCOP
	get_menv_params_for_residue(
		core::pose::Pose const & pose,
		core::conformation::Residue const & rsd,
		core::Size atomno
	) const;


private: // helper methods

	/// @brief Versioning
	core::Size version() const override;


private:

	// Store weight when computing derivatives
	mutable core::Real f_elec_lipidlayer_weight_;

};

} // scoring
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_hh
