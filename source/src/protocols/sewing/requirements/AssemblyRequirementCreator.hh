// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/AssemblyRequirementCreator.hh
/// @brief  Base class for AssemblyRequirementCreators for the AssemblyRequirement load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_sewing_requirements_AssemblyRequirementCreator_hh
#define INCLUDED_protocols_sewing_requirements_AssemblyRequirementCreator_hh

// Unit Headers
#include <protocols/sewing/requirements/AssemblyRequirementCreator.fwd.hh>

// Package Headers
#include <protocols/sewing/requirements/AssemblyRequirement.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace sewing  {
namespace requirements {

/// @brief The Creator class is responsible for creating a particular
/// GlobalRequirement class.
class AssemblyRequirementCreator : public utility::pointer::ReferenceCount
{
public:
	AssemblyRequirementCreator() {}
	virtual ~AssemblyRequirementCreator() {}

	virtual AssemblyRequirementOP create_requirement() const = 0;
	virtual std::string keyname() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const =0;
	//static std::string assembly_requirement_ct_namer( std::string tag_name );

};


} //namespace
} //namespace
} //namespace

#endif
