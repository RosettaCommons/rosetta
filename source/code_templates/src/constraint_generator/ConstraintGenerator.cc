// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

protocols::constraint_generator::ConstraintGeneratorOP
--class--Creator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new --class-- );
}

std::string
--class--Creator::keyname() const
{
	return --class--::class_name();
}

--class--::--class--():
	protocols::constraint_generator::ConstraintGenerator( --class--::class_name() )
{
}

--class--::~--class--()
{}

protocols::constraint_generator::ConstraintGeneratorOP
--class--::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new --class--( *this ) );
}

std::string
--class--::class_name()
{
	return "--class--";
}

void
--class--::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
}

core::scoring::constraints::ConstraintCOPs
--class--::apply( core::pose::Pose const & ) const
{
	return core::scoring::constraints::ConstraintCOPs();
}

void
--class--Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	--class--::provide_xml_schema( xsd );
}

void
--class--::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"DOCUMENTATION STRING FOR CLASS",
		attlist );
}

--end_namespace--
