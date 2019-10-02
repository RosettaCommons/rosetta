// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/EtableEnergy.hh
/// @brief  Etable energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_methods_DistanceChainbreakEnergy_hh
#define INCLUDED_core_scoring_methods_DistanceChainbreakEnergy_hh

// Unit headers
#include <core/energy_methods/DistanceChainbreakEnergy.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//#include <ObjexxFCL/FArray3D.hh>

namespace core {
namespace scoring {
namespace methods {

/// @brief DistanceChainbreakEnergy class iterates across all residues in finalize()
/// and determines the penalty between residues i and i+1 by how far apart
/// their N and C atom are
class DistanceChainbreakEnergy : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:

	DistanceChainbreakEnergy();

	/// clone
	EnergyMethodOP
	clone() const override
	{
		return utility::pointer::make_shared< DistanceChainbreakEnergy >();
	}

	/// called at the end of energy evaluation
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override;
	core::Size version() const override;


};

} // methods
} // scoring
} // core


#endif
