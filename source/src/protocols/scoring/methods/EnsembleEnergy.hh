// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EnsembleEnergy.hh
/// @brief  Function to constrain particular residues to an interface during docking
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_scoring_methods_EnsembleEnergy_hh
#define INCLUDED_protocols_scoring_methods_EnsembleEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace protocols {
namespace scoring {
namespace methods {


class EnsembleEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;
public:


	EnsembleEnergy();


	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	// determine the site_constraint score for all sites constrained to the interface
	virtual
	void
	finalize_total_energy(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const;

	virtual
	core::Size version() const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

};


}
}
}

#endif
