// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsExplicitDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsExplicitDefiner.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_definers/util.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script

using std::string;
using std::endl;
using std::stringstream;
using core::Real;
using core::pose::Pose;
using basic::datacache::DataMap;
using utility::tag::TagCOP;
using utility::vector0;
using basic::Tracer;

namespace protocols {
namespace loops {
namespace loops_definers {


static Tracer TR("protocols.loops.loops_definers.LoopsExplicitDefiner");


LoopsExplicitDefiner::LoopsExplicitDefiner() :
	loop_list_()
{}

LoopsExplicitDefiner::~LoopsExplicitDefiner() {}

LoopsExplicitDefiner::LoopsExplicitDefiner(LoopsExplicitDefiner const & src) : LoopsDefiner(src),
	loop_list_(src.loop_list_)
{}

/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsExplicitDefiner::clone(
) const {
	return LoopsDefinerOP( new LoopsExplicitDefiner(*this) );
}

SerializedLoop
LoopsExplicitDefiner::parse_loop_tag(
	TagCOP const tag,
	string const & loops_name
) {

	SerializedLoop loop;

	if ( tag->hasOption("start") ) {
		loop.start = tag->getOption<Size>("start");
	} else {
		stringstream err_msg;
		err_msg
			<< "Tag " << tag->getName() << " with name "
			<< "'" << loops_name << "' does not have the expected 'start' field." << endl;
		utility_exit_with_message(err_msg.str());
	}

	if ( tag->hasOption("stop") ) {
		loop.stop = tag->getOption<Size>("stop");
	} else {
		stringstream err_msg;
		err_msg
			<< "Tag " << tag->getName() << " with name "
			<< "'" << loops_name << "' does not have the expected 'stop' field." << endl;
		utility_exit_with_message(err_msg.str());
	}

	loop.cut = tag->getOption<Size>("cut", 0);
	loop.skip_rate = tag->getOption<Real>("skip_rate", 0.0);
	loop.extended = tag->getOption<bool>("extended", false);

	return loop;
}


/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsExplicitDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &,
	Pose const &
) {

	if ( !tag->hasOption("name") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, 
			"Unable to create unnamed LoopsDefiner (type: Loops)" );
	}
	string const loops_name(tag->getOption<string>("name"));


	vector0< TagCOP >::const_iterator begin=tag->getTags().begin();
	vector0< TagCOP >::const_iterator end=tag->getTags().end();

	for ( ; begin != end; ++begin ) {
		TagCOP loop_tag= *begin;

		if ( loop_tag->getName() != "loop" ) {
			TR.Error
				<< "Please include only tags with name 'loop' "
				<< "as subtags of a 'Loops' tag" << endl
				<< "Tag with name '" << loop_tag->getName() << "' is invalid" << endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "");
		}

		SerializedLoop loop = parse_loop_tag(loop_tag, loops_name);
		loop_list_.push_back(loop);
	}
}

SerializedLoopList
LoopsExplicitDefiner::apply(
	Pose const &
) {
	return loop_list_;
}

std::string LoopsExplicitDefiner::class_name()
{
	return "Loops";
}

void LoopsExplicitDefiner::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;

	AttributeList loop_attributes;
	loop_attributes
		+ Attr::required_attribute( "start", xsct_non_negative_integer, "The residue index (pose numbering) for the first position in the loop." )
		+ Attr::required_attribute( "stop", xsct_non_negative_integer, "The residue index (pose numbering) for the last position in the loop." )
		+ Attr::attribute_w_default( "cut", xsct_non_negative_integer, "The residue index (pose numbering) for the cutpoint residue (i.e. the cut"
		" residue will be the lower end of the chain break, and the cut+1 residue will be the upper end of the chain break). If a zero value is"
		" given, it will be interpretted as saying that the cutpoint should be chosen using the logic in protocols::loops::Loop::choose_cutpoint", "0" )
		+ Attr::attribute_w_default( "skip_rate", xsct_real, "The probability that you would like the loop will be skipped by the IndependentLoopMover;"
		" a value less than zero means never skip, a value greater than one means always skip", "0.0" )
		+ Attr::attribute_w_default( "extended", xsct_rosetta_bool, "If set to true, then the initial backbone dihedrals for the loop will be overwritten"
		" to produce an extended conformation (i.e. phi=-150 degres, psi=150, omego=180) and the bond angles and lengths will be idealized", "false" );

	XMLSchemaSimpleSubelementList loop_subelements;
	loop_subelements.add_simple_subelement( "loop", loop_attributes, "Define a loop based on start and stop indices" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( class_name() )
		.complex_type_naming_func( & complex_type_name_for_loop_definer )
		.description( "Use this element to define a series of loops in subtags, one loop per subtag" )
		.set_subelements_repeatable( loop_subelements )
		.add_required_name_attribute( "Each Loop element must be given a name for error-reporting purposes, but the name is also used to add the loop collection defined within it to the DataMap by the LoopDefinerLoader" )
		.write_complex_type_to_schema( xsd );

}



} //namespace
} //namespace
} //namespace
