// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loops_definers/util.hh
/// @brief Utility functions useful in LoopDefiner classes.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_util_HH
#define INCLUDED_protocols_loops_loops_definers_util_HH

// Package headers
#include <protocols/loops/Loops.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// basic headers
#include <basic/datacache/DataMap.fwd.hh>

// utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace loops {
namespace loops_definers {

LoopsOP
load_loop_definitions(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	core::pose::Pose const & pose
);

std::string
complex_type_name_for_loop_definer( std::string const & element_name );

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & loop_definer_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);

/// @brief Appends the attributes read by load_loop_definitions
void
attributes_for_load_loop_definitions( utility::tag::AttributeList & attributes );

} //namespace
} // namespace
} // namespace

#endif  // include guard
