// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RMS_Energy.hh
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


#ifndef INCLUDED_core_scoring_methods_RMS_Energy_hh
#define INCLUDED_core_scoring_methods_RMS_Energy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class RMS_Energy : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:


	RMS_Energy();

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
	pose::Pose native_pose_;
	core::Real rms_target_;
	virtual
	core::Size version() const;
};


} // methods
} // scoring
} // core

#endif
