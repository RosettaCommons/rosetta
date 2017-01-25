// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/pose_outputters/pose_outputter_schemas.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_pose_outputters_pose_outputter_schemas_HH
#define INCLUDED_protocols_jd3_pose_outputters_pose_outputter_schemas_HH

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief Define the XML schema definition for a PoseOutputter that contains
/// no subtags but may contain any number of attributes (aka options).
void
pose_outputter_xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & outputter_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);


}  // namespace pose_outputters
}  // namespace jd3
}  // namespace protocols

#endif
