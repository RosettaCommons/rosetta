// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/ReadResfileFromDB.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/task_operations/ReadResfileFromDB.hh>
#include <protocols/task_operations/ReadResfileFromDBCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
//#include <basic/resource_manager/ResourceManager.hh>
//#include <basic/resource_manager/util.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// C++
#include <string>
#include <sstream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace task_operations {

static basic::Tracer TR( "protocols.toolbox.task_operations.ReadResfileFromDB" );

using namespace core::pack::task::operation;
using namespace utility::tag;

using basic::database::get_db_session;
using basic::database::check_statement_sanity;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using core::pose::Pose;
using core::pack::task::parse_resfile_string;
using core::pack::task::ResfileReaderException;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using cppdb::result;
using cppdb::statement;
using std::endl;
using std::string;
using std::stringstream;
using utility::sql_database::sessionOP;
using utility::tag::TagCOP;

TaskOperationOP ReadResfileFromDBCreator::create_task_operation() const {
	return TaskOperationOP( new ReadResfileFromDB );
}

void ReadResfileFromDBCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReadResfileFromDB::provide_xml_schema( xsd );
}

std::string ReadResfileFromDBCreator::keyname() const
{
	return ReadResfileFromDB::keyname();
}

ReadResfileFromDB::ReadResfileFromDB() :
	parent(),
	database_table_("resfiles"),
	db_session_()
{}

ReadResfileFromDB::ReadResfileFromDB(
	utility::sql_database::sessionOP db_session,
	string const & database_table) :
	parent(),
	database_table_(database_table),
	db_session_(std::move(db_session))
{}

ReadResfileFromDB::ReadResfileFromDB(ReadResfileFromDB const & /*src*/) = default;

ReadResfileFromDB::~ReadResfileFromDB() = default;

TaskOperationOP ReadResfileFromDB::clone() const {
	return TaskOperationOP( new ReadResfileFromDB( *this ) );
}

void
ReadResfileFromDB::apply( Pose const & pose, PackerTask & task ) const {

	std::string selection_tag( selection_tag_ );
	if ( selection_tag == "" ) {
		selection_tag = protocols::jd2::current_input_tag();
	}

	stringstream sql_stmt;
	sql_stmt
		<< "SELECT resfile FROM " << database_table_
		<< " WHERE tag='" << selection_tag << "';";
	string sql(sql_stmt.str());
	check_statement_sanity(sql);
	statement select_stmt(safely_prepare_statement(sql, db_session_));
	result res(safely_read_from_database(select_stmt));
	if ( !res.next() ) {
		stringstream error_message;
		error_message
			<< "Unable to locate resfile for job distributor input tag '"
			<< selection_tag << "' in the database." << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str());
	}
	string resfile;
	res >> resfile;
	TR << "Retrieved the following resfile for selection_tag " << selection_tag << "\n";
	TR << resfile;
	TR << std::endl;
	try {
		parse_resfile_string(pose, task, "<resfile from database table '" + database_table_ + "' with tag '" + selection_tag + "'>", resfile);
	} catch( ResfileReaderException const & e ){
		stringstream error_message;
		error_message
			<< "Failed to process resfile stored for input tag '" << selection_tag << "'" << endl
			<< "RESFILE:" << endl
			<< resfile << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str());
	}
}

void
ReadResfileFromDB::db_session(
	utility::sql_database::sessionOP db_session
) {
	db_session_ = db_session;
}

void
ReadResfileFromDB::database_table(string const & database_table) {
	database_table_ = database_table;
}

std::string const &
ReadResfileFromDB::database_table() const {
	return database_table_;
}

void
ReadResfileFromDB::selection_tag( std::string const & setting )
{
	selection_tag_ = setting;
}

std::string const &
ReadResfileFromDB::selection_tag() const
{
	return selection_tag_;
}


void
ReadResfileFromDB::parse_tag( TagCOP tag , DataMap & datamap )
{
	if ( tag->hasOption("database_table") ) {
		database_table_ = tag->getOption<string>("database_table");
	} else if ( tag->hasOption("table") ) {
		database_table_ = tag->getOption<string>("table");
	}

	if ( ! tag->hasOption( "selection_tag" ) ) {
		std::ostringstream oss;
		oss << "ERROR: The ReadResfileFromDB task operation requires a selection_tag attribute to identify which row of the table to pull results from\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}
	selection_tag( tag->getOption< std::string >( "selection_tag" ));
	db_session_ = protocols::rosetta_scripts::parse_database_session(tag, datamap);
}

void ReadResfileFromDB::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "database_table", xs_string , "The table in the database from which to read" )
		+ XMLSchemaAttribute( "table", xs_string , "The table in the database from which to read (same as database_table)" )
		+ XMLSchemaAttribute::required_attribute( "selection_tag", xs_string, "The tag to use to identify"
		" the row in the indicated table that will be read from for the indicated job. In JD3, this"
		" can/should be combined with the script_vars flag so that different jobs can read different"
		" resfiles. This is a marked departure from the JD2 functionality which relied on the global"
		" data representing the currently-running job. That functionality is now removed." );

	protocols::rosetta_scripts::attributes_for_parse_database_session( xsd, attributes );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "This task operation will query a"
		" database for a resfile string. The table in the database should have a column named 'resfile'"
		" and a (key) column named 'tag'. The 'selection_tag' attribute for this task operation should"
		" be used to say which row to query from the database." );
}


} //namespace task_operations
} //namespace protocols
