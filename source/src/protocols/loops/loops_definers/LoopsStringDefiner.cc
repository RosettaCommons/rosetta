// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsStringDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsStringDefiner.hh>
#include <protocols/loops/loops_definers/LoopsStringDefinerCreator.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_definers/util.hh>
#include <protocols/loops/util.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script

using std::string;
using std::stringstream;
using core::Real;
using core::pose::Pose;
using basic::datacache::DataMap;
using utility::tag::TagCOP;
using basic::Tracer;

namespace protocols {
namespace loops {
namespace loops_definers {


static Tracer TR("protocols.loops.loops_definers.LoopsStringDefiner");


LoopsStringDefiner::LoopsStringDefiner(std::string const & in) :
	loop_spec_(in)
{}

LoopsStringDefiner::~LoopsStringDefiner() = default;

LoopsStringDefiner::LoopsStringDefiner(LoopsStringDefiner const & /*src*/) = default;

/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsStringDefiner::clone(
) const {
	return utility::pointer::make_shared< LoopsStringDefiner >(*this);
}

/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsStringDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &
) {

	if ( !tag->hasOption("name") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,
			"Unable to create unnamed LoopsDefiner (type: Loops)" );
	}
	string const loops_name(tag->getOption<string>("name"));

	if ( tag->hasOption("loops") ) {
		loop_spec_ = tag->getOption<string>("loops");
	} else {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Tag with name '"+loops_name+"' does not have the expected 'loops' field.");
	}
}

SerializedLoopList
LoopsStringDefiner::apply(
	Pose const & pose
) {
	return serialized_loops_from_string( loop_spec_, pose );
}

std::string LoopsStringDefiner::class_name()
{
	return "LoopsString";
}

void LoopsStringDefiner::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute( "loops", utility::tag::xs_string , "The format for loops is: Start:End:Cut,Start:End:Cut... RosettaResnum or PDB Numbering accepted.");

	xsd_type_definition_w_attributes( xsd, class_name(), "Define a set of loops by a simple string format.", attributes );
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Creator

LoopsStringDefinerCreator::LoopsStringDefinerCreator() = default;
LoopsStringDefinerCreator::~LoopsStringDefinerCreator() = default;

LoopsDefinerOP
LoopsStringDefinerCreator::create_loops_definer() const {
	return utility::pointer::make_shared< LoopsStringDefiner >();
}

string
LoopsStringDefinerCreator::type_name() const {
	return LoopsStringDefiner::class_name();
}

void LoopsStringDefinerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopsStringDefiner::provide_xml_schema( xsd );
}


} //namespace
} //namespace
} //namespace
