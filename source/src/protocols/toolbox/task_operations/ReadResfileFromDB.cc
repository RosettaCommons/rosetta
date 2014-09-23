// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ReadResfileFromDB.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/toolbox/task_operations/ReadResfileFromDB.hh>
#include <protocols/toolbox/task_operations/ReadResfileFromDBCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

// C++
#include <string>
#include <sstream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

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
using protocols::jd2::JobDistributor;
using std::endl;
using std::string;
using std::stringstream;
using utility::sql_database::sessionOP;
using utility::tag::TagCOP;

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
	db_session_(db_session)
{}

ReadResfileFromDB::ReadResfileFromDB(ReadResfileFromDB const & src) : TaskOperation(src),
		database_table_(src.database_table_),
		db_session_(src.db_session_)
{}

ReadResfileFromDB::~ReadResfileFromDB() {}

TaskOperationOP ReadResfileFromDBCreator::create_task_operation() const {
	return TaskOperationOP( new ReadResfileFromDB );
}

TaskOperationOP ReadResfileFromDB::clone() const {
	return TaskOperationOP( new ReadResfileFromDB( *this ) );
}

void
ReadResfileFromDB::apply( Pose const & pose, PackerTask & task ) const {

	string tag(JobDistributor::get_instance()->current_job()->input_tag());

	stringstream sql_stmt;
	sql_stmt
		<< "SELECT resfile FROM " << database_table_
		<< " WHERE tag='" << tag << "';";
	string sql(sql_stmt.str());
	check_statement_sanity(sql);
	statement select_stmt(safely_prepare_statement(sql, db_session_));
	result res(safely_read_from_database(select_stmt));
	if(!res.next()){
		stringstream error_message;
		error_message
			<< "Unable to locate resfile for job distributor input tag '"
			<< tag << "' in the database." << endl;
		throw utility::excn::EXCN_Msg_Exception(error_message.str());
	}
	string resfile;
	res >> resfile;
	try{
		parse_resfile_string(pose, task, resfile);
	} catch(ResfileReaderException e){
		stringstream error_message;
		error_message
			<< "Failed to process resfile stored for input tag '" << tag << "'" << endl
			<< "RESFILE:" << endl
			<< resfile << endl;
		throw utility::excn::EXCN_Msg_Exception(error_message.str());
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
ReadResfileFromDB::parse_tag( TagCOP tag , DataMap & )
{
	using namespace basic::resource_manager;

	if(tag->hasOption("db")){
		utility_exit_with_message(
			"The 'db' tag has been deprecated. Please use 'database_name' instead.");
	}

	if(tag->hasOption("db_mode")){
		utility_exit_with_message(
			"The 'db_mode' tag has been deprecated. "
			"Please use the 'database_mode' instead.");
	}

	if(tag->hasOption("database_table")){
		database_table_ = tag->getOption<string>("database_table");
	} else if(tag->hasOption("table")){
		database_table_ = tag->getOption<string>("table");
	}

	if(tag->hasOption("resource_description")){
		std::string resource_description = tag->getOption<string>("resource_description");
		if(!ResourceManager::get_instance()->has_resource_with_description(
				resource_description))
		{
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a resource_description of " + resource_description +
					" for ReadResfileFromDB, but the ResourceManager doesn't have a resource with that description" );
		}
		db_session_ = get_resource< utility::sql_database::session >( resource_description );
	}
	else{
		db_session_ = basic::database::parse_database_connection(tag);
	}

}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
