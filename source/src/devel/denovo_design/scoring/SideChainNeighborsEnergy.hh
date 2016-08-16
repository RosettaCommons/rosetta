// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/SideChainNeighborsEnergy.hh
/// @brief  energy term use for scoring predicted burial
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_core_scoring_methods_SideChainNeighborsEnergy_hh
#define INCLUDED_core_scoring_methods_SideChainNeighborsEnergy_hh

#include <devel/denovo_design/scoring/SideChainNeighborsEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class SideChainNeighborsEnergy : public ContextDependentOneBodyEnergy {
public:
	// default constructor
	SideChainNeighborsEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	void
	setup_for_scoring(
		pose::Pose & pose, ScoreFunction const &
	) const;

	void
	finalize_total_energy(
		pose::Pose &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	virtual
	core::Size version() const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const;

private:
	utility::vector1< core::Real > pred_burial_;
};

}
}
}

#endif // INCLUDED_core_scoring_methods_SideChainNeighborsEnergy_HH
