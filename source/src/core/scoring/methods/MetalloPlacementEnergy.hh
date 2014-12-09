// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MetalloPlacementEnergy.hh
/// @brief  Low-res placement score for metal sites
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_methods_MetalloPlacementEnergy_hh
#define INCLUDED_core_scoring_methods_MetalloPlacementEnergy_hh

// Unit Headers
#include <core/scoring/methods/MetalloPlacementEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/SecondaryStructurePotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers


namespace core {
namespace scoring {
namespace methods {

///
class MetalloPlacementEnergy : public WholeStructureEnergy {

	core::Real collision_thresh2_, cb_cb_dis_thresh2_;

public:

	///
	MetalloPlacementEnergy();

	// ///
	// MetalloPlacementEnergy( MetalloPlacementEnergy const & src );
	//

	/// clone
	virtual
	EnergyMethodOP
	clone() const;


	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

virtual
void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {}

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
