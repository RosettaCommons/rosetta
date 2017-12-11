// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loops_definers/util.cc
/// @brief Utility functions useful in LoopDefiner classes.
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

// Project headers
#include <core/pose/Pose.hh>

// basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace loops {
namespace loops_definers {

LoopsOP
load_loop_definitions(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	core::pose::Pose const & pose
) {
	using namespace loops;
	using namespace loops_definers;

	std::string loops_str( tag->getOption< std::string >( "loops" ) );

	if ( data.has("loops_definers", loops_str) ) {
		LoopsDefinerOP loops_definer(
			data.get_ptr< LoopsDefiner >("loops_definers", loops_str));
		return LoopsOP( new Loops(loops_definer->apply(pose)) );
	} else {
		return loops_from_string(loops_str, pose );
	}

	return 0;

}

std::string
complex_type_name_for_loop_definer( std::string const & element_name )
{
	return "loop_definer_" + element_name + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & loop_definer_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	// If the developer has already added a name attribute, do not add another one.
	utility::tag::AttributeList local_attrs( attributes );
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "name", local_attrs ) ) {
		local_attrs + utility::tag::optional_name_attribute();
	}

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_loop_definer )
		.element_name( loop_definer_type )
		.description( description )
		.add_attributes( local_attrs )
		.write_complex_type_to_schema( xsd );
}

//loops( loops::loops_definers::load_loop_definitions(tag, data, pose) )
/// @brief Appends the attributes read by load_loop_definitions

void
attributes_for_load_loop_definitions( utility::tag::AttributeList & attributes )
{
	attributes + utility::tag::XMLSchemaAttribute( "loops", utility::tag::xs_string , "The format for loops is: Start:End:Cut,Start:End:Cut... ");
}

} // namespace
} // namespace
} // namespace
