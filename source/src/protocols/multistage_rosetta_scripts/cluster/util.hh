// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/util.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_util_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_util_HH

#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <string>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

inline std::string
complex_type_name_for_cluster_metric( std::string const & cm_type ){
	return "cm_" + cm_type + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & cm_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);

} // cluster
} // multistage_rosetta_scripts
} // utility

#endif
