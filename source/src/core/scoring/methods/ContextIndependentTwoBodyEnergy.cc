// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextIndependentTwoBodyEnergy.cc
/// @brief  Short-ranged, context-independent, two-body energy class implementation
/// @author Phil Bradley

// Unit Headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

ContextIndependentTwoBodyEnergy::ContextIndependentTwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

ContextIndependentTwoBodyEnergy::~ContextIndependentTwoBodyEnergy() = default;

EnergyMethodType
ContextIndependentTwoBodyEnergy::method_type() const
{
	return ci_2b;
}

}
}
}

