// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/BurialEnergy.hh
/// @brief  energy term use for scoring predicted burial
/// @author James Thompson

#ifndef INCLUDED_core_energy_methods_BurialEnergy_hh
#define INCLUDED_core_energy_methods_BurialEnergy_hh

#include <core/energy_methods/BurialEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class BurialEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy {
public:

	BurialEnergy() : core::scoring::methods::ContextDependentOneBodyEnergy( utility::pointer::make_shared< BurialEnergyCreator >() ) {
		init_from_file();
	}

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	void
	setup_for_scoring(
		pose::Pose & pose, core::scoring::ScoreFunction const &
	) const override;

	void
	finalize_total_energy(
		pose::Pose &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override {}

	core::Size version() const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override;

private:
	void init_from_file();
	utility::vector1< core::Real > pred_burial_;
};

}
}

#endif // INCLUDED_core_energy_methods_BurialEnergy_HH
