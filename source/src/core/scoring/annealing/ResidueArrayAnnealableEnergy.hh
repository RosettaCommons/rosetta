// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/annealing/ResidueArrayAnnealableEnergy.hh
/// @brief  Annealable method interrface for score types evaluated over explicit list of residues.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_scoring_annealing_ResidueArrayAnnealableEnergy_hh
#define INCLUDED_core_scoring_annealing_ResidueArrayAnnealableEnergy_hh

#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace annealing {

class ResidueArrayAnnealableEnergy
{
public:

	/// @brief Constructor.
	///
	ResidueArrayAnnealableEnergy();

	/// @brief Copy constructor.
	///
	ResidueArrayAnnealableEnergy( ResidueArrayAnnealableEnergy const &src );

	/// @brief Destructor.
	///
	virtual ~ResidueArrayAnnealableEnergy();

	/// @brief Calculate the energy given a vector of const-owning pointers to Residue objects.
	/// @details Must be implemented by derived classes.
	virtual core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const = 0;

	/// @brief ResidueArrayAnnealableEnergy objects may optionally cache data within the EnergyMethod prior to a packer run.
	/// This function is defined as doing nothing by default, but can be redefined on a per-EnergyMethod basis to cache whatever
	/// data are necessary.
	/// @details Note that this is strictly unidirectional: data can be cached FROM the Pose or ScoreFunction, not TO the Pose or ScoreFunction.
	/// The method is nonconst to permit caching within the ResidueArrayAnnealableEnergy-derived EnergyMethod.
	virtual void setup_residuearrayannealableenergy_for_packing ( core::pose::Pose const &pose, core::scoring::ScoreFunction const &sfxn);

};

}
}
}

#endif
