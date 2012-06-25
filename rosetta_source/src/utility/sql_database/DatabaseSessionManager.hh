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
/// @author Chris Miles


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

class session : public cppdb::session, public utility::pointer::ReferenceCount {

public:
	void
	set_db_mode(
		DatabaseMode::e const db_mode) { db_mode_ = db_mode; }

	DatabaseMode::e
	get_db_mode() const { return db_mode_; }

private:
	DatabaseMode::e db_mode_;

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
		std::string const & db_name,
		std::string const & pq_schema,
		std::string const & host,
		std::string const & user,
		std::string const & password,
		platform::Size port,
		bool readonly = false,
		bool separate_db_per_mpi_process = false);


	///@brief Acquire a sqlite3 database session
	sessionOP
	get_session_sqlite3(
		std::string const & database,
		bool const readonly=false,
		bool const separate_db_per_mpi_process=false);

	///@brief Acquire a mysql database session
	sessionOP
	get_session_mysql(
		std::string const & database,
		std::string const & host,
		std::string const & user,
		std::string const & password,
		platform::Size port);


	///@brief Acquire a postgres database session
	sessionOP
	get_session_postgres(
		std::string const & database,
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
