// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorRegistrator.hh
/// @brief  Declaration of the template class for registrating ConstraintGeneratorCreators with
///         the ConstriantGeneratorFactory
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_constraint_generator_ConstraintGeneratorRegistrator_hh
#define INCLUDED_protocols_constraint_generator_ConstraintGeneratorRegistrator_hh

// Package headers
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace protocols {
namespace constraint_generator {

/// @brief This templated class will register an instance of an
/// ConstraintGeneratorCreator (class T) with the ConstraintGeneratorFactory.  It will ensure
/// that no ConstraintGenerator creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ConstraintGeneratorRegistrator : public utility::factory::WidgetRegistrator< ConstraintGeneratorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ConstraintGeneratorFactory, T > parent;
public:
	ConstraintGeneratorRegistrator() : parent() {}
};

}
}

#endif
