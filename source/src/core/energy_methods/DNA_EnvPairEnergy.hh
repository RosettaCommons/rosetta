// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/DNA_EnvPairEnergy.hh
/// @brief  dna scoring
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_DNA_EnvPairEnergy_HH
#define INCLUDED_core_scoring_methods_DNA_EnvPairEnergy_HH

// Unit Headers
#include <core/energy_methods/DNA_EnvPairEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/dna/DNA_EnvPairPotential.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

namespace core {
namespace scoring {
namespace methods {

/// @brief  Implementation of env and pair terms for protein-DNA interactions
/// @details  Could be a CI2B, but centroid atom is not currently the nbr atom for dna so intxn threshold tricky

class DNA_EnvPairEnergy : public WholeStructureEnergy {
public:
	typedef WholeStructureEnergy  parent;

public:

	///
	DNA_EnvPairEnergy();


	/// clone
	EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////


	/// @brief  All the work happens here
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & emap
	) const override;


	/// @brief  No graphs required.
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {};

	///
	Size
	version() const override;


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	dna::DNA_EnvPairPotential const & potential_;

};


}
}
}

#endif
