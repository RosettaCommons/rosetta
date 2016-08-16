// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/MPCbetaEnergy.hh
///
/// @brief  Membrane Environemnt CBeta Energy
/// @details One Body Term - Score packing density in the membrane. Scores centroids for within
///    6A and 12A radius. Derived from Membrane base potential and uses mpframework data
///    Last Modified: 4/2/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPCbetaEnergy_hh
#define INCLUDED_core_scoring_membrane_MPCbetaEnergy_hh

// Unit Headers
#include <core/scoring/membrane/MPCBetaEnergy.fwd.hh>

// Package Headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane Environemnt CBeta Energy Term
class MPCbetaEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy {

public:
	typedef ContextDependentOneBodyEnergy  parent;

public: // constructors

	/// @brief Default Constructor for CBeta energy
	MPCbetaEnergy();

	/// @brief Clone Energy Method
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

public: // scoring methods

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &
	) const;

	virtual
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	// Get version
	virtual
	core::Size version() const;

public: // energy function methods

	/// @brief Energy Function for CBeta Term
	core::Real
	compute_mpcbeta_score(
		CenListInfo const & cenlist,
		core::Size const seqpos,
		core::Size const num_tmh
	) const;


private: //data

	// MP Potential Base Instance
	MembraneData const & mpdata_;

};

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPCbetaEnergy_hh
