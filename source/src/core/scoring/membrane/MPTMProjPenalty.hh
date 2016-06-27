// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPTMProjPenalty.hh
///
/// @brief  Membrane Protein TM Proj Penalty
/// @details Whole structure energy - Penalty for unreasonable tm-helix length compared to predicted
///    helix length (from topology) and uses mpframework data
///    Last Modified: 4/3/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTMProjPenalty_hh
#define INCLUDED_core_scoring_membrane_MPTMProjPenalty_hh

// Unit Headers
#include <core/scoring/membrane/MPTMProjPenalty.fwd.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

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

namespace core {
namespace scoring {
namespace membrane {

/// @brief Class Membrane TM proj Penalty
class MPTMProjPenalty : public methods::WholeStructureEnergy {

public: // typedefs

	typedef methods::WholeStructureEnergy  parent;

public: // constructors

	/// @brief Default Constructor
	MPTMProjPenalty();

	/// @brief Clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

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

	/// @brief Compute Penalty for Length of Helix
	core::Real
	compute_tmproj_penalty(
		core::Real start_z_pos,
		core::Real end_z_pos,
		core::Real dist
	) const;

private:

	/// @brief Version
	core::Size version() const { return (core::Size)2.0; }

	// MP Base potential (database)
	MembraneData const & mpdata_;

}; // MPTMProjPenalty

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPTMProjPenalty_hh
