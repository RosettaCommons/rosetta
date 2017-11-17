// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Protocol headers
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/moves/MoverFactory.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/Pose.hh>

// C++ headers
#include <string>
#include <sstream>

using namespace std;

using utility::tag::TagCOP;
using utility::pointer::dynamic_pointer_cast;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using core::pose::Pose;

namespace protocols {
namespace loop_modeling {
namespace utilities {

LoopMoverOP loop_mover_from_tag(
	TagCOP tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose) {

	using protocols::moves::MoverOP;
	using protocols::moves::MoverFactory;

	MoverOP base_mover = MoverFactory::get_instance()->newMover(
		tag, data, filters, movers, pose);
	LoopMoverOP loop_mover = dynamic_pointer_cast< LoopMover > ( base_mover );

	if ( ! loop_mover ) {
		stringstream message;
		message << "<" << tag->getName() << "> is not a loop mover, so ";
		message << "cannot be used in the <LoopModeler> protocol." << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
	}

	return loop_mover;
}

/// @brief return a Loops object from tags
/// @details returns a LoopsOP object from tags "loops_file" or a Loop subtag.  Moved from loop_modeling/LoopMover.cc.
protocols::loops::LoopsOP parse_loops_from_tag( utility::tag::TagCOP tag ) {
	protocols::loops::LoopsOP parsed_loops;

	if ( tag->hasOption("loops_file") ) {
		string loops_file = tag->getOption<string>("loops_file");
		parsed_loops = LoopsOP( new Loops(loops_file) );
	}

	for ( utility::tag::TagCOP subtag : tag->getTags("Loop") ) {
		core::Size start = subtag->getOption<Size>("start");
		core::Size stop = subtag->getOption<Size>("stop");
		core::Size cut = subtag->getOption<Size>("cut", 0);
		core::Real skip_rate = subtag->getOption<Real>("skip_rate", 0.0);
		bool extended = subtag->getOption<bool>("rebuild", false);

		if ( ! parsed_loops ) {
			parsed_loops = LoopsOP( new Loops() );
		}
		parsed_loops->add_loop(start, stop, cut, skip_rate, extended);
	}

	return parsed_loops;  //note may be NULL if no loops were found!
}

///@brief util function for append_subelement_and_attributes_for_parse_loops_from_tag
std::string loop_subelement_namer( std::string const & subelement_name ) { return "loop_subelement_" + subelement_name + "_type"; }
///@brief util function for append_subelement_and_attributes_for_parse_loops_from_tag
//std::string unnamed_loop_ct_namer( std::string const & subelement_name ) { return "unnamed_loop" + subelement_name + "_type"; }

///@brief use with XML schema if you use parse_loops_from_tag
///@details appends the loopfile attribute to the AttributeList for the top level thing (the Mover itself), and appends the schema for Loops subtags to the provided (usually empty, but not necessarily) subelements list
void append_subelement_and_attributes_for_parse_loops_from_tag (
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::XMLSchemaSimpleSubelementList & subelements,
	utility::tag::AttributeList & attributes )
{
	using namespace utility::tag;

	//First, append the attribute for the optional loops_file element
	attributes + XMLSchemaAttribute( "loops_file", xs_string, "path to loops file" );

	//This generates schema for the Loop subtag / subelement

	//Some of this is copied from LoopsExplicitDefiner's provide_xml_schema
	//AttributeList for Loop element
	AttributeList loop_attributes;
	loop_attributes
		+ XMLSchemaAttribute::required_attribute( "start", xsct_non_negative_integer, "The residue index (pose numbering) for the first position in the loop." )
		+ XMLSchemaAttribute::required_attribute( "stop", xsct_non_negative_integer, "The residue index (pose numbering) for the last position in the loop." )
		+ XMLSchemaAttribute::attribute_w_default( "cut", xsct_non_negative_integer, "The residue index (pose numbering) for the cutpoint residue (i.e. the cut"
		" residue will be the lower end of the chain break, and the cut+1 residue will be the upper end of the chain break). If a zero value is"
		" given, it will be interpretted as saying that the cutpoint should be chosen using the logic in protocols::loops::Loop::choose_cutpoint", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "skip_rate", xsct_real, "The probability that you would like the loop will be skipped (such as by the IndependentLoopMover);"
		" a value less than zero means never skip, a value greater than one means always skip - NOT RESPECTED BY ALL PROTOCOLS", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rebuild", xsct_rosetta_bool, "If set to true, then the initial backbone dihedrals for the loop will be overwritten"
		" to produce an extended conformation (i.e. phi=-150 degres, psi=150, omego=180) and the bond angles and lengths will be idealized - NOT RESPECTED BY AL PROTOCOLS", "false" );

	//Generate type for Loop element
	XMLSchemaComplexTypeGenerator unnamed_loop;
	std::string const loop("Loop");
	unnamed_loop.element_name( loop )
		.description( "Use this element to define a series of loops in subtags, one loop per subtag" )
		.complex_type_naming_func( & loop_subelement_namer )
		.add_attributes(loop_attributes)
		.write_complex_type_to_schema( xsd );  //this last call adds it to the global schema (safely ignoring repeated adds), but does not modify the schema of the class calling this function

	//this call modifies the schema of the class calling this function
	subelements.add_already_defined_subelement( loop, & loop_subelement_namer );

	return;
}

/// @brief For XML schema: Appends the attributes read by set_task_factory_from_tag, which simply calls "parse_task_operations"
void
attributes_for_set_task_factory_from_tag( utility::tag::AttributeList & attributes )
{
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attributes );
}

/// @brief For XML scheema: Appends the attributes read by set_scorefxn_from_tag, which simply calls  "parse_score_function"
void
attributes_for_set_scorefxn_from_tag( utility::tag::AttributeList & attributes )
{
	protocols::rosetta_scripts::attributes_for_parse_score_function( attributes );
}

}
}
}
