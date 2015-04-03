// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextDependentTwoBodyEnergy.cc
/// @brief  Score function class
/// @author Phil Bradley

// Unit Headers
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers


namespace core {
namespace scoring {
namespace methods {

ContextDependentTwoBodyEnergy::ContextDependentTwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

ContextDependentTwoBodyEnergy::~ContextDependentTwoBodyEnergy() {}


EnergyMethodType
ContextDependentTwoBodyEnergy::method_type() const
{
	return cd_2b;
}

}
}
}

