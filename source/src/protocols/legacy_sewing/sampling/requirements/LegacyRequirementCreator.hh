// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/RequirementCreator.hh
/// @brief  Base class for RequirementCreators for the Requirement load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementCreator_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementCreator_hh

// Unit Headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementCreator.fwd.hh>

// Package Headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// c++ headers
#include <string>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

/// @brief The Creator class is responsible for creating a particular
/// LegacyGlobalRequirement class.
class LegacyGlobalRequirementCreator : public utility::pointer::ReferenceCount
{
public:
	LegacyGlobalRequirementCreator() {}
	virtual ~LegacyGlobalRequirementCreator() {}

	virtual LegacyGlobalRequirementOP create_requirement() const = 0;
	virtual std::string type_name() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

class LegacyIntraSegmentRequirementCreator : public utility::pointer::ReferenceCount
{
public:
	LegacyIntraSegmentRequirementCreator() {}
	virtual ~LegacyIntraSegmentRequirementCreator() {}

	virtual LegacyIntraSegmentRequirementOP create_requirement() const = 0;
	virtual std::string type_name() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

} //namespace
} //namespace
} //namespace
} //namespace

#endif
