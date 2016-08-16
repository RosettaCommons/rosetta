// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MembraneEnvPenalties.hh
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


#ifndef INCLUDED_core_scoring_methods_MembraneEnvPenalties_hh
#define INCLUDED_core_scoring_methods_MembraneEnvPenalties_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class MembraneEnvPenalties : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:


	MembraneEnvPenalties();

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

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:
	// const-ref to scoring database
	MembranePotential const & potential_;
	virtual
	core::Size version() const;
};


} // methods
} // scoring
} // core

#endif
