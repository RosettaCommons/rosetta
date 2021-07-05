// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/CenPairEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_energy_methods_CenPairEnergy_hh
#define INCLUDED_core_energy_methods_CenPairEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class CenPairEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;
public:


	CenPairEnergy();


	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &, // pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// totals
	) const override;

	/// This method *should* admit to defining intraresidue energies
	bool
	defines_intrares_energy( core::scoring::EnergyMap const & ) const override { return false; }

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override {}

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	core::scoring::EnvPairPotential const & potential_;
	core::Size version() const override;

};


}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
