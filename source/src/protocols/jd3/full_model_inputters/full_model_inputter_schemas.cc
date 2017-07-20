// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/full_model_inputters/full_model_inputter_schemas.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit header
#include <protocols/jd3/full_model_inputters/full_model_inputter_schemas.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputterFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

void
full_model_inputter_xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & inputter_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & FullModelInputterFactory::complex_type_name_for_full_model_inputter )
		.element_name( inputter_type )
		.description( description )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}


}  // namespace full_model_inputters
}  // namespace jd3
}  // namespace protocols
