// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/ContactOrderEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author James Thompson


#ifndef INCLUDED_core_energy_methods_ContactOrderEnergy_hh
#define INCLUDED_core_energy_methods_ContactOrderEnergy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class ContactOrderEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;
public:

	ContactOrderEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	core::Real
	calculate_contact_order( pose::Pose const & pose ) const;
	core::Size version() const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:


}; // ContactOrderEnergy


} // scoring
} // core

#endif // INCLUDED_core_energy_methods_ContactOrderEnergy_HH
