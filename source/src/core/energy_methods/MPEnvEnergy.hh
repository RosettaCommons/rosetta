// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/MPEnvEnergy.hh
///
/// @brief  Membrane Environemnt Energy
/// @details One Body Term - score residue interaction with specific hydrophobic layer
///    derived from Membrane base potential and uses mpframework data
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPEnvEnergy_hh
#define INCLUDED_core_scoring_membrane_MPEnvEnergy_hh

// Unit Headers
#include <core/energy_methods/MPEnvEnergy.fwd.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>

#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/membrane/MembraneData.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/AA.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Mmebrane Environemnt Energy Term
class MPEnvEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy {

public: // typedefs

	typedef core::scoring::methods::ContextDependentOneBodyEnergy parent;

public: // methods

	/// @brief Constructor
	MPEnvEnergy();

	/// @brief Clone Method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief Setup Menv for Scoring by updating neighbor count and compute cen env
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	/// @brief Setup for Derivatives - Calls Setup for Scoring
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const override;

	/// @brief Compute Membrane Environemnt - Residue Energy
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const override;


	/// @brief Finalize Whole Structure Energy from One Body Energies
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &
	) const override;

	/// @brief No Required Context Graphs in MP Terms
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}

public: // energy method functions

	/// @brief Full Evaluate Membrane Env Method (given z_position and const menv score)
	core::Real
	compute_mpenv_score(
		CenListInfo const & cenlist,
		core::chemical::AA const & aa,
		core::Real const z_position,
		core::Size const seqpos
	) const;
private:

	// MP Env and Cen Rotamer Statistics
	MembraneData const & mpdata_;

	// Versioning
	core::Size version() const override;

}; // MPEnvEnergy

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPEnvEnergy_hh
