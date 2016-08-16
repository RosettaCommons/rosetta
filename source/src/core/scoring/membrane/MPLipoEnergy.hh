// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/MPLipoEnergy.hh
///
/// @brief  Membrane Lipophibicity Term
/// @details Whole Structure Energy - Evaluate structure based on derived
///    lipophobicities from input in lips file.
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPLipoEnergy_hh
#define INCLUDED_core_scoring_membrane_MPLipoEnergy_hh

// Unit Headers
#include <core/scoring/membrane/MPLipoEnergy.fwd.hh>

// Package Headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/membrane/MembraneData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane Lipophilicity term
class MPLipoEnergy : public core::scoring::methods::WholeStructureEnergy  {

public:

	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:

	/// Energy Method Creators /////////////////

	/// @brief Defalt Constructor
	MPLipoEnergy();

	/// @brief Clone Method
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/// Scoring Methods ////////////////////////

	/// @brief Score WHole Structure Energy
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	/// @brief Setup for Scoring
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const {}

private:

	// const-ref to scoring database
	MembraneData const & mpdata_;

	// Include lips in scoring
	bool include_lips_;

	// versioning
	virtual
	core::Size version() const;
};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPLipoEnergy_hh
