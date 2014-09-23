// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/DatabaseSessionManager.hh
/// @author Matthew O'Meara
/// @author Sam Deluca
/// @author Tim Jacobs

#ifndef INCLUDED_utility_sql_database_DatabaseSessionManager_HH
#define INCLUDED_utility_sql_database_DatabaseSessionManager_HH

// Unit Headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/sql_database/types.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <platform/types.hh>

// Boost Headers
#include <boost/scoped_ptr.hpp>

// C++ Headers
#include <string>

// External
#include <cppdb/frontend.h>

namespace utility {
namespace sql_database {

//TODO alexford
// Need to refactor this class to move data members into correct members of 
// cppdb::session. Specifically, should generate a connection specific information
// class to store connection metadata.


//@brief cppdb::session subclass containing connection data, only pass via reference.
class session : public cppdb::session, public utility::pointer::ReferenceCount {

public:
	session() :
		cppdb::session(),
		db_mode_(DatabaseMode::sqlite3),
		db_name_(""),
		pq_schema_(""),
		db_partition_(-1),
		cur_transaction_(/* 0 */),
		chunk_size_(0),
		transaction_counter_(0)
	{
	}

	~session();

	void
	set_db_mode(
		DatabaseMode::e const db_mode) { db_mode_ = db_mode; }

	DatabaseMode::e
	get_db_mode() const { return db_mode_; }

	void
	set_transaction_mode(
		TransactionMode::e const transaction_mode) { transaction_mode_ = transaction_mode; }

	TransactionMode::e
	get_transaction_mode() const { return transaction_mode_; }

	void
	set_chunk_size(
		platform::Size const chunk_size) { chunk_size_ = chunk_size; }

	platform::Size
	get_chunk_size() { return chunk_size_; }

	void
	set_db_name(
		std::string const & db_name) { db_name_ = db_name; }

	std::string const &
	get_db_name() const { return db_name_; }

	void
	set_pq_schema(
		std::string const & pq_schema) { pq_schema_ = pq_schema; }

	std::string const &
	get_pq_schema() const { return pq_schema_; }

	platform::SSize get_db_partition() { return db_partition_; }
	void set_db_partition(platform::SSize db_partition) { db_partition_ = db_partition; }
	bool is_db_partitioned() { return db_partition_ > 0; }

	///@brief indicate that a transaction block has begun.
	void begin_transaction();

	///@brief indicate that a transaction block has completed. NOTE:
	///When in chunk transaction mode, this may not actually write to
	///the database.
	void
	commit_transaction();

	///@brief force a transaction to be committed. This should only
	///be used when writing data that is required by other processes,
	///such as protocol and batch ids.
	void
	force_commit_transaction();

private:
	DatabaseMode::e db_mode_;
	std::string db_name_;
	std::string pq_schema_;

	platform::SSize db_partition_;

	transactionOP cur_transaction_;
	TransactionMode::e transaction_mode_;
	platform::Size chunk_size_;
	platform::Size transaction_counter_;

};

//@brief cppdb::transaction subclass for smart pointers.
class transaction : public cppdb::transaction, public utility::pointer::ReferenceCount {

public:

	transaction(session & session) :
		cppdb::transaction(session)
	{}

};

class DatabaseSessionManager {

protected:

	// Private constructor to make it singleton managed
	DatabaseSessionManager();
	DatabaseSessionManager( const DatabaseSessionManager & src );

public:

	// Warning this is not called because of the singleton pattern
	virtual ~DatabaseSessionManager();

	///@brief return singleton instance of session manager
	static
	DatabaseSessionManager *
	get_instance();

	///@brief Acquire a database session
	sessionOP
	get_db_session(
		DatabaseMode::e db_mode,
		TransactionMode::e transaction_mode,
		platform::Size chunk_size,
		std::string const & db_name,
		std::string const & pq_schema,
		std::string const & host,
		std::string const & user,
		std::string const & password,
		platform::Size port,
		bool readonly = false,
		platform::SSize db_partition = -1);


	///@brief Acquire a sqlite3 database session
	sessionOP
	get_session_sqlite3(
		std::string const & database,
		TransactionMode::e transaction_mode=TransactionMode::standard,
		platform::Size chunk_size=0,
		bool const readonly=false,
		platform::SSize db_partition=-1);

	///@brief Acquire a mysql database session
	sessionOP
	get_session_mysql(
		std::string const & database,
		TransactionMode::e transaction_mode,
		platform::Size chunk_size,
		std::string const & host,
		std::string const & user,
		std::string const & password,
		platform::Size port);


	///@brief Acquire a postgres database session
	sessionOP
	get_session_postgres(
		std::string const & database,
		TransactionMode::e transaction_mode,
		platform::Size chunk_size,
		std::string const & pq_schema,
		std::string const & host,
		std::string const & user,
		std::string const & password,
		platform::Size port);

private:

	void
	set_postgres_schema_search_path(
		sessionOP db_session,
		utility::vector1< std::string > const & schema_search_path
	);



#ifndef MULTITHREADED
	static boost::scoped_ptr< DatabaseSessionManager > instance_;
#endif
};



}
}



#endif
