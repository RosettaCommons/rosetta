// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorCreator.hh
/// @brief  Class for instantiating a particular ConstraintGenerator
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_constraint_generator_ConstraintGeneratorCreator_HH
#define INCLUDED_protocols_constraint_generator_ConstraintGeneratorCreator_HH

// Package headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// C++ headers
#include <string>

namespace protocols {
namespace constraint_generator {

class ConstraintGeneratorCreator : public utility::pointer::ReferenceCount {
public:
	/// @brief Instantiate a particular ConstraintGenerator
	virtual ConstraintGeneratorOP
	create_constraint_generator() const = 0;

	/// @brief Return a string that will be used to instantiate the particular ConstraintGenerator
	virtual std::string
	keyname() const = 0;
	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const =0;
};


} //namespace constraint_generator
} //namespace protocols


#endif
