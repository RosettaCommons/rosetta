// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPPairEnergy.hh
///
/// @brief  Membrane Residue Pair Energy Term
/// @details Two Body Term - score residue-residue interactions in the membrane. Derived from Membrane
///    base potential and uses mpframework data
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPPairEnergy_hh
#define INCLUDED_core_scoring_membrane_MPPairEnergy_hh

// Unit Headers
#include <core/scoring/membrane/MPPairEnergy.hh>

// Package headers
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/membrane/MembraneData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/chemical/AA.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// Rebecca, our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::scoring;

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane Environemtn Residue Pair Energy Term
class MPPairEnergy : public core::scoring::methods::ContextDependentTwoBodyEnergy {

public:
	typedef ContextDependentTwoBodyEnergy parent;

public: // constructors

	/// @brief Default Constructor
	MPPairEnergy();

	/// @brief Clone Method
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

public: // scoring methods

	/// @brief Setup for Scoring - compute cen env and update neighbors
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/// @brief Compute residue pair energy in th emembrane
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	/// @brief Fianlize whole structure energy
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & ) const { return false; }

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	/// @brief Define Atomic Interaction Cutoff == 6A
	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	/// @brief Versioning
	virtual
	core::Size version() const;

public: // energy methods

	/// @brief Compute Reisdue-Residue Pair energy from z_position, aa, and Cb-Cb distance
	core::Real
	compute_mpair_score(
		core::Real const z_position1,
		core::Real const z_position2,
		chemical::AA const & aa1,
		chemical::AA const & aa2,
		core::Real const cendist
	) const;

private: // data

	// MP Potential Base Instance
	MembraneData const & mpdata_;

	// User option
	bool no_interpolate_mpair_;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPPairEnergy_hh
