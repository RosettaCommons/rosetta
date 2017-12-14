// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/parser/JumpSelectorLoader.hh>
#include <protocols/parser/JumpSelectorLoaderCreator.hh>

// Project headers
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/select/jump_selector/JumpSelectorFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.JumpSelectorLoader" );

JumpSelectorLoader::JumpSelectorLoader() = default;
JumpSelectorLoader::~JumpSelectorLoader() = default;

void JumpSelectorLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	using namespace utility::tag;
	using core::select::jump_selector::JumpSelectorOP;
	using TagCOPs = utility::vector0<TagCOP>;

	TagCOPs const & selector_tags( tag->getTags() );
	for ( core::Size ii = 0; ii < selector_tags.size(); ++ii ) {
		TagCOP ii_tag = selector_tags[ ii ];
		JumpSelectorOP selector = core::select::jump_selector::JumpSelectorFactory::get_instance()->new_jump_selector(
			ii_tag->getName(),
			ii_tag,
			datamap
		);

		// if "name" is specified, add it to the data map under that name. Otherwise use the type name.
		bool const data_add_status = datamap.add( "JumpSelector" ,
			ii_tag->getOption( "name", ii_tag->getName() ),
			selector );
		if ( !data_add_status ) {
			utility_exit_with_message( "JumpSelector '" + ii_tag->getName() + "' already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
	TR.flush();
}

std::string
JumpSelectorLoader::loader_name() { return "JUMP_SELECTORS"; }

std::string
JumpSelectorLoader::jump_selector_loader_ct_namer( std::string const & element_name )
{
	return "jump_selector_loader_" + element_name + "_type";
}

void JumpSelectorLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace core::select::jump_selector;

	JumpSelectorFactory::get_instance()->define_jump_selector_xml_schema( xsd );

	XMLSchemaSimpleSubelementList js_loader_subelements;
	js_loader_subelements.add_group_subelement( & JumpSelectorFactory::jump_selector_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator js_ct;
	js_ct.element_name( loader_name() ).complex_type_naming_func( & jump_selector_loader_ct_namer )
		.description( "JumpSelectors may be defined as subelements of the " + loader_name() + " element, and then will be placed into the DataMap"
		" for later retrieval by Movers, Filters, MoveMapFactories or anything else that might use a JumpSelector. All immediate subelements should have the 'name' attribute"
		" as that is how they will be identified in the DataMap. Subelements of the immediate subelements will not be loaded into the data map and do"
		" not need to be given a name; e.g. if an immediate subelemement is an 'And' JumpSelector, and this 'And' JumpSelector has a 'Chain' subelement,"
		" then the 'And' subelement must be given a name, but the 'Chain' subelement probably should not be given a name. Why not? The 'Chain' subelement"
		" will not end up as a member of the DataMap in any case, and if a name is given to it, then a user of the XML file may think it is reasonable"
		" to use the name for that Chain selector elsewhere in their script -- but their script would fail if they do so." )
		.set_subelements_repeatable( js_loader_subelements )
		.write_complex_type_to_schema( xsd );

}


DataLoaderOP
JumpSelectorLoaderCreator::create_loader() const { return DataLoaderOP( new JumpSelectorLoader ); }

std::string
JumpSelectorLoaderCreator::keyname() const { return JumpSelectorLoader::loader_name(); }

JumpSelectorLoaderCreator::DerivedNameFunction
JumpSelectorLoaderCreator::schema_ct_naming_function() const
{
	return & JumpSelectorLoader::jump_selector_loader_ct_namer;
}

void
JumpSelectorLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	JumpSelectorLoader::provide_xml_schema( xsd );
}



} //namespace parser
} //namespace protocols
