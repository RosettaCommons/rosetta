// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file	 utility/sql_database/DatabaseSessionManager.cc
/// @author Matthew O'Meara
/// @author Sam Deluca
/// @author Chris Miles

#ifdef USEMPI
#include <mpi.h>
#endif

// Unit Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/util.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/assert.hh>
#include <platform/types.hh>

// Boost Headers
#ifdef MULTITHREADED
#include <boost/thread/tss.hpp>
#else
#include <boost/scoped_ptr.hpp>
#endif

// C++ Headers
#include <string>
#include <sstream>

namespace utility {
namespace sql_database {

using utility::file::FileName;
using std::string;
using std::stringstream;
using cppdb::cppdb_error;
using platform::Size;
using cppdb::statement;

#ifdef MULTITHREADED
boost::thread_specific_pointer< DatabaseSessionManager > DatabaseSessionManager::instance_;
#else
boost::scoped_ptr< DatabaseSessionManager > DatabaseSessionManager::instance_;
#endif

DatabaseSessionManager *
DatabaseSessionManager::get_instance(){
	if( instance_.get() == 0 ){
		instance_.reset( new DatabaseSessionManager() );
	}
	return instance_.get();
}

DatabaseSessionManager::DatabaseSessionManager() {}

DatabaseSessionManager::DatabaseSessionManager(
const DatabaseSessionManager &
) {}

DatabaseSessionManager::~DatabaseSessionManager() {}

sessionOP
DatabaseSessionManager::get_db_session(
	DatabaseMode::e db_mode,
	string const & db_name,
	string const & pq_schema,
	string const & host,
	string const & user,
	string const & password,
	Size port,
	bool readonly,
	bool separate_db_per_mpi_process
){

	switch(db_mode){
	case DatabaseMode::sqlite3:
		return get_session_sqlite3(db_name, readonly, separate_db_per_mpi_process);
	case DatabaseMode::mysql:
		return get_session_mysql(db_name, host, user, password, port);
	case DatabaseMode::postgres:
		return get_session_postgres(
			db_name, pq_schema, host, user, password, port);
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_mode) + "'");
	}
}


///details@ For SQLite3 database, the separate_db_per_mpi_process
/// appends "_<mpi_rank>" to the end of the database filename This is
/// useful when writing to an sqlite database not through the job
/// distributor where locking causes problems
sessionOP
DatabaseSessionManager::get_session_sqlite3(
	string const & database,
	bool const readonly /* = false */,
	bool const MPI_ONLY( separate_db_per_mpi_process ) /* = false */
){
	sessionOP s(new session());
	s->set_db_mode(DatabaseMode::sqlite3);

#ifdef USEMPI
	string use_database;
	if(separate_db_per_mpi_process){
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		stringstream buf; buf << FileName(database).name() << "_" << mpi_rank;
		use_database = buf.str();
	} else {
		use_database = FileName(database).name();
	}
#else
	string const use_database(FileName(database).name());
#endif


	try {
		if(readonly){
			s->open("sqlite3:mode=readonly;db="+use_database);
		} else {
			s->open("sqlite3:db="+use_database);
		}
	} catch (cppdb_error & e){
		std::stringstream error_msg;
		error_msg
			<< "Failed to open sqlite3 database file '" << database << "'"
			<< (readonly ? " in readonly mode:" : ":") << std::endl
			<< "\t" << e.what();
		utility_exit_with_message(error_msg.str());
	}
	return s;
}

sessionOP
DatabaseSessionManager::get_session_mysql(
	string const & MYSQL_ONLY(database),
	string const & MYSQL_ONLY(host),
	string const & MYSQL_ONLY(user),
	string const & MYSQL_ONLY(password),
  Size MYSQL_ONLY(port)
){

#ifndef USEMYSQL
	utility_exit_with_message(
		"If you want to use a mysql database, build with extras=mysql");
	return 0;
#else

	sessionOP s(new session());
	s->set_db_mode(DatabaseMode::mysql);

	stringstream connection_string;
	connection_string
		<< "mysql:host=" << host << ";"
		<< "user=" << user << ";"
		<< "password=" << password << ";"
		<< "database=" << database << ";"
		<< "port=" << port << ";"
		<< "opt_reconnect=1";

	try {
		s->open(connection_string.str());
	} catch (cppdb_error & e){
		std::stringstream error_msg;
		error_msg
			<< "Failed to open mysql database:"
			<< "\thost='" << host << "'" << std::endl
			<< "\tuser='" << user << "'" << std::endl
			<< "\tpassword='**********'" << std::endl
			<< "\tport='" << port << "'" << std::endl
			<< "\tdatabase='" << database << "'" << std::endl
			<< std::endl
			<< "\t" << e.what();
		utility_exit_with_message(error_msg.str());
	}
	return s;

#endif
}

sessionOP
DatabaseSessionManager::get_session_postgres(
	string const & POSTGRES_ONLY(database),
	string const & POSTGRES_ONLY(pq_schema),
	string const & POSTGRES_ONLY(host),
	string const & POSTGRES_ONLY(user),
	string const & POSTGRES_ONLY(password),
  Size POSTGRES_ONLY(port)
){

#ifndef USEPOSTGRES
	utility_exit_with_message(
		"If you want to use a postgres database, build with extras=postgres");
	return 0;
#else

	sessionOP s(new session());
	s->set_db_mode(DatabaseMode::postgres);

	stringstream connection_string;
	connection_string
		<< "postgresql:host=" << host << ";"
		<< "user=" << user << ";"
		<< "password=" << password << ";"
		<< "port=" << port << ";"
		<< "dbname=" << database;

	try {
		s->open(connection_string.str());
	} catch (cppdb_error & e){
		std::stringstream error_msg;
		error_msg
			<< "Failed to open postgres database:"
			<< "\thost='" << host << "'" << endl
			<< "\tuser='" << user << "'" << endl
			<< "\tpassword='**********'" << endl
			<< "\tport='" << port << "'" << endl
			<< "\tdatabase='" << database << "'" << endl
			<< endl
			<< "\t" << e.what();
		utility_exit_with_message(error_msg.str());
	}

	vector1<string> schema_search_path;
	schema_search_path.push_back(pq_schema);
	set_postgres_schema_search_path(s, schema_search_path);

	return s;

#endif
}

///@detail postgres does not allow queries between databases, instead
/// it allows tables to be created in different namespaces called
/// "schemas". By specifing the search path, statements will be
/// executed in a specified namespace. Note setting the search path
/// only affects this session.
///
/// For example, to use the schema UBQdesign_stage1_r456644_120323,
/// set set the search path with the vector ["UBQdesign_stage1_r456644_120323"].
/// See: http://www.postgresql.org/docs/8.1/static/ddl-schemas.html
void
DatabaseSessionManager::set_postgres_schema_search_path(
	sessionOP db_session,
	vector1< string > const & schema_search_path
) {
	stringstream stmt_str;
	stmt_str << "SET search_path TO";
	for(Size i=1; i <= schema_search_path.size(); ++i){
		stmt_str << (i==1 ? "" : ", ") << schema_search_path[i];
	}
	stmt_str << ";";
	statement stmt(db_session->prepare(stmt_str.str()));
	stmt.exec();
}



} // namespace
} // namespace

