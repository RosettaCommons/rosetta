// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/WholeStructureEnergy.cc
/// @brief  Implementation of the EnergyMethod that applies to a whole structure and not on a residue-by-residue basis
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

WholeStructureEnergy::WholeStructureEnergy(
	EnergyMethodCreatorOP creator
) : parent( creator )
{}

EnergyMethodType
WholeStructureEnergy::method_type() const
{
	return ws;
}

} // methods
} // scoring
} // core
