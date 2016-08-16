// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentOneBodyEnergy.cc
/// @brief
/// @author Phil Bradley

// Unit headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

ContextIndependentOneBodyEnergy::ContextIndependentOneBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}


EnergyMethodType
ContextIndependentOneBodyEnergy::method_type() const
{
	return ci_1b;
}

} // methods
} // scoring
} // core
