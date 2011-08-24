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
#include <core/pose/Pose.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

// C++
#include <string>
#include <sstream>

namespace protocols {
namespace toolbox {
namespace task_operations {

using basic::database::get_db_session;
using core::pose::Pose;
using core::pack::task::parse_resfile_string;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using cppdb::result;
using cppdb::statement;
using protocols::jd2::JobDistributor;
using std::endl;
using std::string;
using std::stringstream;
using utility::sql_database::sessionOP;
using utility::tag::TagPtr;

ReadResfileFromDB::ReadResfileFromDB() :
	parent(),
	database_filename_("resfiles.db3"),
	database_mode_("sqlite3"),
	database_table_("resfiles")
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if(option[inout::database_filename].user()){
		database_filename_ = option[inout::database_filename].value();
		database_mode_ = option[inout::database_mode].value();
	}
}

ReadResfileFromDB::ReadResfileFromDB(
	string const & database_filename,
	string const & database_mode,
	string const & database_table) :
	parent(),
	database_filename_(database_filename),
	database_mode_(database_mode),
	database_table_(database_table)
{}

ReadResfileFromDB::ReadResfileFromDB(
	ReadResfileFromDB const & src) :
	database_filename_(src.database_filename_),
	database_mode_(src.database_mode_),
	database_table_(src.database_table_)
{}

ReadResfileFromDB::~ReadResfileFromDB() {}

TaskOperationOP ReadResfileFromDBCreator::create_task_operation() const {
	return new ReadResfileFromDB;
}

TaskOperationOP ReadResfileFromDB::clone() const {
	return new ReadResfileFromDB( *this );
}

void
ReadResfileFromDB::apply( Pose const &, PackerTask & task ) const {

	string tag(JobDistributor::get_instance()->current_job()->input_tag());

	sessionOP db_session(
		get_db_session(database_filename_, database_mode_, true));

	stringstream sql_stmt;
	sql_stmt
		<< "SELECT resfile FROM " << database_table_
		<< " WHERE tag='" << tag << "';";
	result res = (*db_session) << sql_stmt.str();
	if(!res.next()){
		stringstream error_message;
		error_message
			<< "Unable to locate resfile for job distributor input tag '"
			<< tag << "' in database '" << database_filename_ << "'." << endl;
		utility_exit_with_message(error_message.str());
	}
	string resfile;
	res >> resfile;
	parse_resfile_string(task, resfile);
}

void
ReadResfileFromDB::database_filename(string const & database_filename) {
	database_filename_ = database_filename;
}

std::string const &
ReadResfileFromDB::database_filename() const {
	return database_filename_;
}

void
ReadResfileFromDB::database_mode(string const & database_mode) {
	database_mode_ = database_mode;
}

std::string const &
ReadResfileFromDB::database_mode() const {
	return database_mode_;
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
ReadResfileFromDB::parse_tag( TagPtr tag )
{
	if(tag->hasOption("db")){
		database_filename_ = tag->getOption<string>("db");
	}
	if(tag->hasOption("db_mode")){
		database_mode_ = tag->getOption<string>("db_mode");
	}

	if(tag->hasOption("table")){
		database_table_ = tag->getOption<string>("table");
	}
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
