// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelectorCreator.hh
/// @brief  Class for instantiating a particular JumpSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_JumpSelectorCreator_HH
#define INCLUDED_core_select_jump_selector_JumpSelectorCreator_HH

// Package headers
#include <core/select/jump_selector/JumpSelector.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace select {
namespace jump_selector {

class JumpSelectorCreator : public utility::pointer::ReferenceCount {
public:
	/// @brief Instantiate a particular JumpSelector
	virtual JumpSelectorOP create_jump_selector() const = 0;

	/// @brief Return a string that will be used to instantiate the particular JumpSelector
	/// from an XML file -- the name for the tag. E.g. "Neighborhood" for the NeighborhoodJumpSelector
	virtual std::string keyname() const = 0;

	/// @brief Define the structure of the XML file for the JumpSelector that this
	/// %JumpSelectorCreator instantiates using the XML Schema language.
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const = 0;
};


} //namespace jump_selector
} //namespace select
} //namespace core


#endif
