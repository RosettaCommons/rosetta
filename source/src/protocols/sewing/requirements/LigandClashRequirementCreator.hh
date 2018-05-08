// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/requirements/LigandClashRequirementCreator.hh
/// @brief  Creates LigandClashRequirements for the AssemblyRequirementFactory
/// @author Minnie Langlois minnie@email.unc.edu

#ifndef INCLUDED_protocols_sewing_requirements_LigandClashRequirementCreator_hh
#define INCLUDED_protocols_sewing_requirements_LigandClashRequirementCreator_hh

// Unit Headers
#include <protocols/sewing/requirements/LigandClashRequirementCreator.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>


namespace protocols {
namespace sewing  {
namespace requirements {

/// @brief The Creator class is responsible for creating a particular
/// Requirement class.
class LigandClashRequirementCreator : public AssemblyRequirementCreator
{
public:
	LigandClashRequirementCreator() {}
	virtual ~LigandClashRequirementCreator() {}

	virtual AssemblyRequirementOP create_requirement() const override;
	virtual std::string keyname() const override;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};


} //namespace
} //namespace
} //namespace

#endif
