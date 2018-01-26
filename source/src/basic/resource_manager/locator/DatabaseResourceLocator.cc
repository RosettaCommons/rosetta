// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/DatabaseResourceLocater.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/locator/DatabaseResourceLocator.hh>
#include <basic/resource_manager/locator/DatabaseResourceLocatorCreator.hh>

// Package headers
#include <basic/resource_manager/locator/locator_schemas.hh>

//project headers
#include <basic/Tracer.hh>
#include <basic/resource_manager/locator/StringResourceStream.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/database/sql_utils.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//External Headers
#include <cppdb/frontend.h>

//C++ headers
#include <istream>
#include <utility>

namespace basic {
namespace resource_manager {
namespace locator {

using utility::tag::TagCOP;
using std::string;
using std::stringstream;
using std::endl;
using std::istream;
using basic::Tracer;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using cppdb::statement;
using cppdb::result;
using utility::sql_database::session;
using utility::sql_database::sessionOP;

static Tracer TR("basic.resource_manager.locator.DatabaseResourceLocator");


///// DatabaseResourceLocatorCreator /////
DatabaseResourceLocatorCreator::DatabaseResourceLocatorCreator() = default;

DatabaseResourceLocatorCreator::~DatabaseResourceLocatorCreator() = default;

ResourceLocatorOP
DatabaseResourceLocatorCreator::create_resource_locator() const {
	return ResourceLocatorOP( new DatabaseResourceLocator );
}

string
DatabaseResourceLocatorCreator::locator_type() const {
	return DatabaseResourceLocator::classname();
}

void
DatabaseResourceLocatorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	return DatabaseResourceLocator::provide_xml_schema( xsd );
}


///// DatabaseResourceLocator /////

DatabaseResourceLocator::DatabaseResourceLocator() : basic::resource_manager::ResourceLocator() {}

DatabaseResourceLocator::DatabaseResourceLocator(
	std::string const & database_session_resource_tag,
	std::string const & sql_command
) : basic::resource_manager::ResourceLocator(),
	database_session_resource_tag_(database_session_resource_tag),
	sql_command_(sql_command)
{}

DatabaseResourceLocator::DatabaseResourceLocator(
	DatabaseResourceLocator const & src
) : basic::resource_manager::ResourceLocator(),
	database_session_resource_tag_(src.database_session_resource_tag_),
	sql_command_(src.sql_command_),
	column_separator_(src.column_separator_)
{}


void
DatabaseResourceLocator::show(
	std::ostream & out
) const {
	out
		<< "DatabaseResourceLocator:" << endl
		<< "\tdatabase_session_resource_tag: '" << database_session_resource_tag_ << "'" << endl
		<< "\tsql_command: '" << sql_command_ << "'" << endl
		<< "\tcolumn_seprator: '" << column_separator_ << "'" << endl;
}

std::string
DatabaseResourceLocator::type() const {
	return classname();
}


DatabaseResourceLocator::~DatabaseResourceLocator() = default;

/// @brief Create a ResourceStream object from the given resource
/// source, so that its stream can be passed to the ResourceLoader
ResourceStreamOP
DatabaseResourceLocator::locate_resource_stream(
	string const & input_id
) const
{

	statement select_stmt(safely_prepare_statement(sql_command_, db_session_ ));
	select_stmt.bind(1, input_id);
	result res(safely_read_from_database(select_stmt));

	if ( !res.next() ) {
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << input_id << "' returned no rows." << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, err_msg.str());
	}


	if ( res.cols() == 0 || res.cols() == -1 ) {
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << input_id << "' returned '" << res.cols() << "' columns." << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, err_msg.str());
	}

	stringstream concatenated_result;
	for ( Size col = 1, ncols = res.cols(); col <= ncols; ++col ) {
		std::string col_val;
		res >> col_val;
		concatenated_result << col_val;
		if ( col < ncols ) {
			concatenated_result << column_separator_;
		}
	}

	if ( res.next() ) {
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << input_id << "' returned more than on row." << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, err_msg.str());
	}


	StringResourceStreamOP string_resource_stream( new StringResourceStream(concatenated_result.str()) );

	return string_resource_stream;
}

void
DatabaseResourceLocator::parse_my_tag(
	TagCOP tag
)
{
	db_session_ = basic::database::parse_database_connection( tag );

	if ( tag->hasOption("sql_command") ) {
		sql_command_ = tag->getOption<string>("sql_command");
		basic::database::check_statement_sanity(sql_command_);
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,
			"The DatabaseResourceLocator requires a 'sql_command' attribute that is an SQL SELECT statement with one parameter '?' and as a key returns a result set with a single column and and a single row containing the data.");
	}

	column_separator_ = tag->getOption<string>("column_separator", "\n");

}

void
DatabaseResourceLocator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
)
{
	using namespace utility::tag;
	AttributeList attrs;
	basic::database::attributes_for_parse_database_connection( attrs, xsd );
	attrs
		+ XMLSchemaAttribute::required_attribute( "sql_command", xs_string, "The command that will"
		" be used to query the database (along with the 'input_id' of the resource) for the string that"
		" will be used to construct a resource. This command should have a single question mark that"
		" the input_id will be substituted for" )
		+ XMLSchemaAttribute( "column_separator", xs_string, "When multiple columns are returned by the query, the 'column_separator' string can be used to separate each column in the final string that is used to construct the Resource" );

	xsd_type_definition_w_attributes( xsd, classname(), "The database resource locator will open a"
		" session with the indicated database and then use the given sql_command to retrieve data from"
		" the database that will then be used to construct a resource", attrs );
}

std::string
DatabaseResourceLocator::classname()
{
	return "DatabaseResourceLocator";
}


} // namespace locator
} // namespace resource_manager
} // namespace basic
