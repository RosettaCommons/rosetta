// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/annealing/ResidueArrayAnnealableEnergy.cc
/// @brief  Annealable method interrface for score types evaluated over explicit list of residues.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace annealing {

/// @brief Constructor.
///
ResidueArrayAnnealableEnergy::ResidueArrayAnnealableEnergy() {}

/// @brief Copy constructor.
///
ResidueArrayAnnealableEnergy::ResidueArrayAnnealableEnergy( ResidueArrayAnnealableEnergy const &/*src*/ )
//TODO -- copy private member variables here.
{}


/// @brief Destructor.
///
ResidueArrayAnnealableEnergy::~ResidueArrayAnnealableEnergy() {}

/// @brief ResidueArrayAnnealableEnergy objects may optionally cache data within the EnergyMethod prior to a packer run.
/// This function is defined as doing nothing by default, but can be redefined on a per-EnergyMethod basis to cache whatever
/// data are necessary.
void
ResidueArrayAnnealableEnergy::setup_residuearrayannealableenergy_for_packing (
	core::pose::Pose const &/*pose*/,
	core::scoring::ScoreFunction const &/*sfxn*/
) {} //Does nothing unless implemented by a derived class.

}
}
}
