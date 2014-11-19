// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPNonHelixPenalty.hh
///
///	@brief		Membrane Protein Non helix in Mmebrane Penalty
///	@details	Whole structure energy - penalty for helices not in the membrane?
///				and uses mpframework data
///				Last Modified: 3/31/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPNonHelixPenalty_hh
#define INCLUDED_core_scoring_membrane_MPNonHelixPenalty_hh

// Unit Headers
#include <core/scoring/membrane/MPNonHelixPenalty.fwd.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/membrane/Span.fwd.hh>
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
	
/// @brief Class Membrane Non Helix Penalty
class MPNonHelixPenalty : public methods::ContextDependentOneBodyEnergy {
	
public: // typedefs
	typedef ContextDependentOneBodyEnergy  parent;
	
public: // constructors
	
	/// @brief Default Constructor
	MPNonHelixPenalty();
	
	/// @brief Clone
	virtual
	EnergyMethodOP
	clone() const;
	
	
	/// @brief Set MP nonhelix Penalty for Scoring
	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const {}
	
	/// @brief Setup MP nonhelix for derivatives
	virtual
	void
	setup_for_derivatives( pose::Pose &, ScoreFunction const & ) const {}
	
	/// @brief Compute nonhelix penalty per-residue
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
	
public: // compute penalty
	
	/// @brief Compute Non Helix in Membrane Score (per-residue)
	core::Real
	compute_nonhelix_penalty( bool tmregion, char secstruc, core::Real z_position ) const;
	
private:
	
	/// @brief Version
	core::Size version() const { return core::Size(2.0); }
	
	// MP Base potential (database)
	MembraneData const & mpdata_;
	
}; // MPNonHelixPenalty
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPNonHelixPenalty_hh
