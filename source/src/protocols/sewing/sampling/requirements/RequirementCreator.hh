// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/feature/RequirementCreator.hh
/// @brief  Base class for RequirementCreators for the Requirement load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_sewing_sampling_requirements_RequirementCreator_hh
#define INCLUDED_protocols_sewing_sampling_requirements_RequirementCreator_hh

// Unit Headers
#include <protocols/sewing/sampling/requirements/RequirementCreator.fwd.hh>

// Package Headers
#include <protocols/sewing/sampling/requirements/GlobalRequirement.fwd.hh>
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

/// @brief The Creator class is responsible for creating a particular
/// GlobalRequirement class.
class GlobalRequirementCreator : public utility::pointer::ReferenceCount
{
public:
	GlobalRequirementCreator() {}
	virtual ~GlobalRequirementCreator() {}

	virtual GlobalRequirementOP create_requirement() const = 0;
	virtual std::string type_name() const = 0;
};

class IntraSegmentRequirementCreator : public utility::pointer::ReferenceCount
{
public:
	IntraSegmentRequirementCreator() {}
	virtual ~IntraSegmentRequirementCreator() {}

	virtual IntraSegmentRequirementOP create_requirement() const = 0;
	virtual std::string type_name() const = 0;
};

} //namespace
} //namespace
} //namespace
} //namespace

#endif
