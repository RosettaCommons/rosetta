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
#include <protocols/parser/DatabaseSessionLoader.hh>
#include <protocols/parser/DatabaseSessionLoaderCreator.hh>

// Basic headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External headers
#include <cppdb/frontend.h>

namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.DatabaseSessionLoader" );

DatabaseSessionLoader::DatabaseSessionLoader() = default;
DatabaseSessionLoader::~DatabaseSessionLoader() = default;

void DatabaseSessionLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;
	using TagCOPs = utility::vector0<TagCOP>;

	TagCOPs const & session_tags( tag->getTags() );

	for ( TagCOP session_tag : session_tags ) {

		utility::sql_database::sessionOP db_session = basic::database::parse_database_connection( session_tag );
		runtime_assert( session_tag->hasOption( "name" ) );
		std::string session_name = session_tag->getOption< std::string >( "name" );
		bool const data_add_status = data.add( "db_sessions" , session_name, db_session );

		if ( !data_add_status ) {
			utility_exit_with_message( "DatabaseSession " + session_name + " already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
}

std::string
DatabaseSessionLoader::loader_name() { return "DATABASE_SESSIONS"; }

std::string
DatabaseSessionLoader::database_session_loader_ct_namer( std::string const & element_name )
{
	return "database_session_loader_" + element_name + "_type";
}

void
DatabaseSessionLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// dbsession attributes
	AttributeList db_session_attributes;
	basic::database::attributes_for_parse_database_connection( db_session_attributes, xsd );
	db_session_attributes + required_name_attribute( "The name given to the database session that Movers and Filters will use to retrieve it from the DataMap" );

	//  complex type
	XMLSchemaSimpleSubelementList loader_subelements;
	loader_subelements.add_simple_subelement( "DatabaseSession", db_session_attributes,
		"Each DatabaseSession object declared below can be accessed from the DataMap by other scriptable elements" );
	XMLSchemaComplexTypeGenerator loader_ct_gen;
	loader_ct_gen
		.element_name( loader_name() )
		.description( "Create DatabaseSession objects to store in the DataMap that can be used by other scriptable elements to read from- or write to a database" )
		.complex_type_naming_func( & database_session_loader_ct_namer )
		.set_subelements_repeatable( loader_subelements )
		.write_complex_type_to_schema( xsd );
}

DataLoaderOP
DatabaseSessionLoaderCreator::create_loader() const { return DataLoaderOP( new DatabaseSessionLoader ); }

std::string
DatabaseSessionLoaderCreator::keyname() const { return DatabaseSessionLoader::loader_name(); }

DatabaseSessionLoaderCreator::DerivedNameFunction
DatabaseSessionLoaderCreator::schema_ct_naming_function() const
{
	return & DatabaseSessionLoader::database_session_loader_ct_namer;
}

void
DatabaseSessionLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DatabaseSessionLoader::provide_xml_schema( xsd );
}


} //namespace parser
} //namespace protocols
