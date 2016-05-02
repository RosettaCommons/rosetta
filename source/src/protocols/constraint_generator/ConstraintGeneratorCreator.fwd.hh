// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorCreator.fwd.hh
/// @brief  Forward declaration of a class that instantiates a particular ConstraintGenerator
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_protocols_constraint_generator_ConstraintGeneratorCreator_FWD_HH
#define INCLUDED_protocols_constraint_generator_ConstraintGeneratorCreator_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace constraint_generator {

class ConstraintGeneratorCreator;

typedef utility::pointer::shared_ptr< ConstraintGeneratorCreator > ConstraintGeneratorCreatorOP;
typedef utility::pointer::shared_ptr< ConstraintGeneratorCreator const > ConstraintGeneratorCreatorCOP;

} //namespace constraint_generator
} //namespace protocols

#endif
