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

#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
	protocols::filters::Filter( "--class--" )
{

}

--class--::~--class--()
{}

void
--class--::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap &
) {

}

protocols::filters::FilterOP
--class--::clone() const
{
	return utility::pointer::make_shared< --class-- >( *this );
}

/// @brief This filter is unpublished.  It returns --name-- as its author.
void
--class--::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"--class--", basic::citation_manager::CitedModuleType::Filter,
		"--name--",
		"TODO: institution",
		"--email--",
		"Wrote the --class--."
		)
	);
}


protocols::filters::FilterOP
--class--::fresh_instance() const
{
	return utility::pointer::make_shared< --class-- >();
}

bool
--class--::apply( core::pose::Pose const & ) const
{
	/*
		NOTE: If you are implementing a filter that computes something, it is STRONGLY
		recommended that you implement this as a simple metric instead.  This allows the
		value's use in more contexts.  For filtering, a simple wrapper around a simple
		metric (which runs the simple metric and returns true or false based on the value
		returned) may be implemented.
	*/
	return true;
}

core::Real
--class--::report_sm( core::pose::Pose const & ) const
{
	return -99999.9;
}

void
--class--::report( std::ostream &, core::pose::Pose const & ) const
{

}

std::string --class--::name() const {
	return class_name();
}

std::string --class--::class_name() {
	return "--class--";
}

void --class--::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"--brief--",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
--class--Creator::create_filter() const
{
	return utility::pointer::make_shared< --class-- >( );
}

std::string
--class--Creator::keyname() const
{
	return --class--::class_name();
}

void --class--Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	--class--::provide_xml_schema( xsd );
}

--end_namespace--
