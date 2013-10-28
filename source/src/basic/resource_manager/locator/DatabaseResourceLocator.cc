// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/DatabaseResourceLocater.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/locator/DatabaseResourceLocator.hh>
#include <basic/resource_manager/locator/DatabaseResourceLocatorCreator.hh>

//project headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/resource_manager/locator/StringResourceStream.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/database/sql_utils.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//External Headers
#include <cppdb/frontend.h>

//C++ headers
#include <istream>
#include <string>

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

static Tracer TR("basic.resource_manager.locator.DatabaseResourceLocator");


///// DatabaseResourceLocatorCreator /////
DatabaseResourceLocatorCreator::DatabaseResourceLocatorCreator() {}

DatabaseResourceLocatorCreator::~DatabaseResourceLocatorCreator() {}

ResourceLocatorOP
DatabaseResourceLocatorCreator::create_resource_locator() const {
	return new DatabaseResourceLocator;
}

string
DatabaseResourceLocatorCreator::locator_type() const {
	return "DatabaseResourceLocator";
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
	return "DatabaseResourceLocator";
}


DatabaseResourceLocator::~DatabaseResourceLocator() {}

/// @brief Create a ResourceStream object from the given resource
/// source, so that its stream can be passed to the ResourceLoader
ResourceStreamOP
DatabaseResourceLocator::locate_resource_stream(
	string const & locator_id
) const {

	if(!ResourceManager::get_instance()->has_resource(database_session_resource_tag_)){
		stringstream err_msg;
		err_msg
			<< "Attempting to locate the data for key='" << locator_id << "' "
			<< "in the DatabaseResourceLocator with tag '" << locator_tag() << "'. "
			<< "However, a database session could not be found "
			<< "because the resource tag given for the database session, "
			<< "'" << database_session_resource_tag_ << "' "
			<< "does not exist in the ResourceManager.";
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}

	session * db_session(dynamic_cast< session * > (
			ResourceManager::get_instance()->find_resource(database_session_resource_tag_)()));

	if(!db_session){
		stringstream err_msg;
		err_msg
			<< "Attempting to locate the data for key='" << locator_id << "' "
			<< "in the DatabaseResourceLocator with tag '" << locator_tag() << "'. "
			<< "However, a database session could not be found "
			<< "because the resource given for the database session tag, "
			<< "'" << database_session_resource_tag_ << "' "
			<< "could is not a DatabaseResourceLocator.";
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}

	statement select_stmt(safely_prepare_statement(sql_command_, db_session));
	select_stmt.bind(1, locator_id);
	result res(safely_read_from_database(select_stmt));

	if(!res.next()) {
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << locator_id << "' returned no rows." << endl;
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}


	if(res.cols() == 0 || res.cols() == -1){
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << locator_id << "' returned '" << res.cols() << "' columns." << endl;
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}

	stringstream concatenated_result;
	for(Size col = 1, ncols = res.cols(); col <= ncols; ++col){
		std::string col_val;
		res >> col_val;
		concatenated_result << col_val;
		if(col < ncols){
			concatenated_result << column_separator_;
		}
	}

	if(res.next()){
		stringstream err_msg;
		err_msg
			<< "In the DatabaseLocator with tag '" << locator_tag() << "', the query:" << endl
			<< sql_command_ << endl
			<< "with parameter '?' <- '" << locator_id << "' returned more than on row." << endl;
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}


	StringResourceStreamOP string_resource_stream(new StringResourceStream(concatenated_result.str()));

	return string_resource_stream;
}

void
DatabaseResourceLocator::parse_my_tag(
	TagCOP tag
) {
	if(tag->hasOption("database_session_tag")){
		database_session_resource_tag_ = tag->getOption<string>("database_session_tag");
	} else {
		throw utility::excn::EXCN_Msg_Exception
			( "The DatabaseResourceLocator requires a 'database_session_tag' that corresponds with a 'DatabaseSession' resource defined in a <Resources/> block.");
	}

	if(tag->hasOption("sql_command")){
		sql_command_ = tag->getOption<string>("sql_command");
		basic::database::check_statement_sanity(sql_command_);
	} else {
		throw utility::excn::EXCN_Msg_Exception (
			"The DatabaseResourceLocator requires a 'sql_command' tag that is an SQL SELECT statement with one parameter '?' and as a key returns a result set with a single column and and a single row containing the data.");
	}

	column_separator_ = tag->getOption<string>("column_separator", "\n");

}

} // namespace locator
} // namespace resource_manager
} // namespace basic
