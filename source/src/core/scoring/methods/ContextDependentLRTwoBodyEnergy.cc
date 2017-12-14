// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextDependentLRTwoBodyEnergy.cc
/// @brief  Context-Dependent, Long-Range, Two-Body Energy class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


// Project Headers

// ObjexxFCL Headers

namespace core {
namespace scoring {
namespace methods {

ContextDependentLRTwoBodyEnergy::ContextDependentLRTwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

ContextDependentLRTwoBodyEnergy::~ContextDependentLRTwoBodyEnergy() = default;

}
}
}

