// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentLRTwoBodyEnergy.cc
/// @brief  Context-Independent, Long-Range, Two-Body Energy class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

ContextIndependentLRTwoBodyEnergy::ContextIndependentLRTwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

ContextIndependentLRTwoBodyEnergy::~ContextIndependentLRTwoBodyEnergy() {}

EnergyMethodType
ContextIndependentLRTwoBodyEnergy::method_type() const
{
	return ci_lr_2b;
}

}
}
}

