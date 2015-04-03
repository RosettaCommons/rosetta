// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 utility/sql_database/DatabaseSessionManager.cc
/// @author Matthew O'Meara
/// @author Sam Deluca
/// @author Chris Miles

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
using std::endl;
using std::stringstream;
using cppdb::cppdb_error;
using platform::Size;
using platform::SSize;
using cppdb::statement;

#ifdef MULTITHREADED
boost::thread_specific_pointer< DatabaseSessionManager > DatabaseSessionManager::instance_;
#else
boost::scoped_ptr< DatabaseSessionManager > DatabaseSessionManager::instance_;
#endif

session::~session(){
	force_commit_transaction();
}

void
session::begin_transaction(){
	if(!cur_transaction_){
		switch(transaction_mode_){
			case(TransactionMode::none):
				//do nothing
				break;
			case(TransactionMode::standard):
				cur_transaction_ = transactionOP( new transaction(*this) );
				break;
			case(TransactionMode::chunk):
				cur_transaction_ = transactionOP( new transaction(*this) );
				break;
			default:
				utility_exit_with_message(
					"Unrecognized transaction mode: '" +
					name_from_transaction_mode(transaction_mode_) + "'");
		}
	}
}

void
session::commit_transaction(){
	if(cur_transaction_){
		switch(transaction_mode_){
			case(TransactionMode::none):
				utility_exit_with_message("You have specified a transaction mode of none, but have an open transaction. Please file a bug");
				break;
			case(TransactionMode::standard):
				cur_transaction_->commit();
				cur_transaction_.reset();
				break;
			case(TransactionMode::chunk):
				if(transaction_counter_==chunk_size_){
					cur_transaction_->commit();
					transaction_counter_=0;
					cur_transaction_.reset();
				}
				else{
					++transaction_counter_;
				}
				break;
			default:
				utility_exit_with_message(
					"Unrecognized transaction mode: '" +
					name_from_transaction_mode(transaction_mode_) + "'");
		}
	}
}

void
session::force_commit_transaction(){
	if(cur_transaction_){
		switch(transaction_mode_){
			case(TransactionMode::none):
				utility_exit_with_message("You have specified a transaction mode of none, but have an open transaction. Please file a bug");
				break;
			case(TransactionMode::standard):
				cur_transaction_->commit();
				cur_transaction_.reset();
				break;
			case(TransactionMode::chunk):
				cur_transaction_->commit();
				cur_transaction_.reset();
				break;
			default:
				utility_exit_with_message(
					"Unrecognized transaction mode: '" +
					name_from_transaction_mode(transaction_mode_) + "'");
		}
	}
}

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
	TransactionMode::e transaction_mode,
	Size chunk_size,
	string const & db_name,
	string const & pq_schema,
	string const & host,
	string const & user,
	string const & password,
	Size port,
	bool readonly,
	SSize db_partition)
{

	switch(db_mode){
	case DatabaseMode::sqlite3:
		return get_session_sqlite3(
							db_name,
							transaction_mode,
							chunk_size,
							readonly,
							db_partition);
	case DatabaseMode::mysql:
		return get_session_mysql(
							db_name,
							transaction_mode,
							chunk_size,
							host,
							user,
							password,
							port);
	case DatabaseMode::postgres:
		return get_session_postgres(
							db_name,
							transaction_mode,
							chunk_size,
							pq_schema,
							host,
							user,
							password,
							port);
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_mode) + "'");
	}
	return 0;
}


///details@ For SQLite3 database, the db_partition is app
/// appends "_<mpi_rank>" to the end of the database filename This is
/// useful when writing to an sqlite database not through the job
/// distributor where locking causes problems
sessionOP
DatabaseSessionManager::get_session_sqlite3(
	string const & database,
	TransactionMode::e transaction_mode,
	Size chunk_size,
	bool const readonly /* = false */,
	SSize const db_partition
){
	sessionOP s( new session() );
	s->set_db_mode(DatabaseMode::sqlite3);
	s->set_transaction_mode(transaction_mode);
	s->set_chunk_size(chunk_size);
	s->set_db_name(database);
	s->set_db_partition(db_partition);

	stringstream connection_string;

	connection_string << "sqlite3:";

	if(readonly)
	{
		connection_string << "mode=readonly;";
	}

	connection_string << "db=" << FileName(database).name();

	if(s->is_db_partitioned())
	{
		connection_string << "_" << s->get_db_partition();
	}

	try {
		s->open(connection_string.str());
	} catch (cppdb_error & e){
		std::stringstream error_msg;
		error_msg
			<< "Failed to open sqlite3 database name '" << database << "'"
			<< " with connection string: " << connection_string.str() << std::endl
			<< "\t" << e.what();
		utility_exit_with_message(error_msg.str());
	}
	return s;
}

sessionOP
DatabaseSessionManager::get_session_mysql(
	string const & MYSQL_ONLY(database),
	TransactionMode::e MYSQL_ONLY(transaction_mode),
	Size MYSQL_ONLY(chunk_size),
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
	s->set_transaction_mode(transaction_mode);
	s->set_chunk_size(chunk_size);
	s->set_db_name(database);

	stringstream connection_string;
	connection_string
		<< "mysql:host=" << host << ";"
		<< "user=" << user << ";"
		<< "password=" << password << ";"
		<< "database=" << database << ";"
		<< "port=" << port << ";"
		<< "opt_reconnect=1";

	platform::Size retry_count = 0;
	platform::Size max_retry = 10;
	//Occasionally a connection will fail to an SQL database due to a busy server,
	//random communications fluke, or something else.  If this happens, try a few more times
	//before giving up entirely.
	while(retry_count < max_retry)
	{
		retry_count++;
		try {
			s->open(connection_string.str());
			break;
		} catch (cppdb_error & e){

			if(retry_count == max_retry)
			{
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
			}else
			{
				sleep(1);
			}
		}
	}


	return s;

#endif
}

sessionOP
DatabaseSessionManager::get_session_postgres(
	string const & POSTGRES_ONLY(database),
	TransactionMode::e POSTGRES_ONLY(transaction_mode),
	Size POSTGRES_ONLY(chunk_size),
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
	s->set_transaction_mode(transaction_mode);
	s->set_chunk_size(chunk_size);
	s->set_db_name(database);
	s->set_pq_schema(pq_schema);

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
			<< "Failed to open postgres database:" << endl
			<< "\thost='" << host << "'" << endl
			<< "\tuser='" << user << "'" << endl
			<< "\tpassword='**********'" << endl
			<< "\tport='" << port << "'" << endl
			<< "\tdatabase='" << database << "'" << endl
			<< "\tpq_schema='" << pq_schema << "'" << endl
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

/// @detail postgres does not allow queries between databases, instead
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
	stmt_str << "SET search_path TO ";
	for(Size i=1; i <= schema_search_path.size(); ++i){
		stmt_str << (i==1 ? "" : ", ") << schema_search_path[i];
	}
	stmt_str << ";";
	statement stmt(db_session->prepare(stmt_str.str()));
	stmt.exec();
}


} // namespace
} // namespace

