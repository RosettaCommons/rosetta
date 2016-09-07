// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/HydrogenBondConstraintGeneratorCreator.hh
/// @brief  Base class for ConstraintCreators for the Constraint
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGeneratorCreator_HH
#define INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGeneratorCreator_HH

// Project Headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

namespace protocols {
namespace constraint_generator {

class HydrogenBondConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	protocols::constraint_generator::ConstraintGeneratorOP create_constraint_generator() const override;
	std::string keyname() const override;

};

} //namespace constraint_generator
} //namespace protocols

#endif
