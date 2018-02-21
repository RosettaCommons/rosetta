// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/util.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/multistage_rosetta_scripts/cluster/util.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & cm_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
) {
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_cluster_metric )
		.element_name( cm_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}


} // cluster
} // multistage_rosetta_scripts
} // utility
