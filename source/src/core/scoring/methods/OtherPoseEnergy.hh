// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/OtherPoseEnergy.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_methods_OtherPoseEnergy_HH
#define INCLUDED_core_scoring_methods_OtherPoseEnergy_HH

#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace scoring {
namespace methods {


class OtherPoseEnergy : public WholeStructureEnergy {

public:
	typedef WholeStructureEnergy  parent;


	OtherPoseEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	virtual
	core::Size version() const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:


};


} //methods
} //scoring
} //core

#endif
