// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacySegmentLengthRequirementCreator.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacySegmentLengthRequirementCreator_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacySegmentLengthRequirementCreator_hh

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementCreator.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.fwd.hh>

//Utility headers
#include <string>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

class LegacySegmentLengthRequirementCreator : public LegacyIntraSegmentRequirementCreator {
	virtual LegacyIntraSegmentRequirementOP create_requirement() const;
	virtual std::string type_name() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};


} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
