// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPTerminiPenalty.hh
///
/// @brief  Membrane Protein Termini Penalty
/// @details Whole structure energy - penalty for residues on the wrong side of the membrane?
///    nd uses mpframework data
///    Last Modified: 3/31/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTerminiPenalty_hh
#define INCLUDED_core_scoring_membrane_MPTerminiPenalty_hh

// Unit Headers
#include <core/scoring/membrane/MPTerminiPenalty.fwd.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {

/// @brief Class Membrane Termini Penalty
class MPTerminiPenalty : public methods::ContextDependentOneBodyEnergy {

public: // typedefs
	typedef ContextDependentOneBodyEnergy  parent;

public: // constructors

	/// @brief Default Constructor
	MPTerminiPenalty();

	/// @brief Clone
	virtual
	EnergyMethodOP
	clone() const;

	/// @brief Set MP Termini Penalty for Scoring
	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const {}

	/// @brief Setup MP termini for derivatives
	virtual
	void
	setup_for_derivatives( pose::Pose &, ScoreFunction const & ) const {}

	/// @brief Compute termini penalty per-residue
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Finalize total energy method (for whole structure
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const {}

public: // penalty method

	/// @brief Compute Termini Penalty
	core::Real
	compute_termini_penalty( core::Real z_position ) const;

private:

	/// @brief Version
	core::Size version() const { return core::Size(2.0); }

	// MP Base potential (database)
	MembraneData const & mpdata_;

}; // MPTerminiPenalty

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPTerminiPenalty_hh
