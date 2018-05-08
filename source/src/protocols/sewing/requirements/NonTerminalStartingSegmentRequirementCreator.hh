// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/requirements/NonTerminalStartingSegmentRequirementCreator.hh
/// @brief  Creates NonTerminalStartingSegmentRequirements for the AssemblyRequirementFactory
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirementCreator_hh
#define INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirementCreator_hh

// Unit Headers
#include <protocols/sewing/requirements/NonTerminalStartingSegmentRequirementCreator.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>

// Package Headers
#include <protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.fwd.hh>

namespace protocols {
namespace sewing  {
namespace requirements {

/// @brief The Creator class is responsible for creating a particular
/// Requirement class.
class NonTerminalStartingSegmentRequirementCreator : public AssemblyRequirementCreator
{
public:
	NonTerminalStartingSegmentRequirementCreator() {}
	virtual ~NonTerminalStartingSegmentRequirementCreator() {}

	virtual AssemblyRequirementOP create_requirement() const;
	virtual std::string
	keyname() const;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;

};


} //namespace
} //namespace
} //namespace

#endif
