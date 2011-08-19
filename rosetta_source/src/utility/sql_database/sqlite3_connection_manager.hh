// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/sqlite3_connection_manager.hh
/// @brief  Easy Access to sqlite3 database; manage connections
/// @author Matthew O'Meara



#ifndef INCLUDE_utility_sql_database_sqlite3_connection_manager_HH
#define INCLUDE_utility_sql_database_sqlite3_connection_manager_HH

#ifdef DB_SQLITE3
// this interface requires the DB_SQLITE3 compilation flag
// Enable with $./scons.py extras=sqlite

// Unit Headers
#include <utility/sql_database/sqlite3_connection_manager.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// External Headers
#include <sqlite3.h>

// C++ Headers
#include <list>
#include <map>
#include <ostream>

namespace utility {
namespace sql_database {


class Sqlite3ConnectionManager : public utility::pointer::ReferenceCount {

public: // typedefs

	///@brief datastructure class to manage connections for each database
	class ConnectionList : public utility::pointer::ReferenceCount {
	public:
		ConnectionList();

		ConnectionList(
			sqlite3 * & db_connection
		);

		virtual ~ConnectionList();

		void
		free_all();

		void
		push_front(
			sqlite3 * & db_connection
		);

		void
		extract_front(
			sqlite3 * & db_connection
		);

		bool
		empty();

		void
		show(
			std::ostream & out
		) const ;

		friend
		std::ostream &
		operator<< (
			std::ostream & out,
			const ConnectionList & connection_list
		);



	private:
		std::list< sqlite3 * > connection_list_;
	};

	typedef pointer::owning_ptr< ConnectionList > ConnectionListOP;

	typedef std::map< std::string , ConnectionListOP > ConnectionLists;



private:


	// Private constructor to implement singleton pattern
	Sqlite3ConnectionManager();

	virtual ~Sqlite3ConnectionManager();

public: // Public Interface

	///@brief return singleton instance of connection manager
	static
	Sqlite3ConnectionManagerOP
	get_instance();

	///@brief request connection to specified database
	///return SQLITE3 result code:
	///    http://www.sqlite.org/c3ref/c_abort.html
	int
	get_connection(
		std::string const & db_fname,
		sqlite3 * & db_connection
	);


	///@brief when done with the connection, free it so others may use it.
	void
	free_connection(
		std::string const db_fname,
		sqlite3 * & db_connection
	);


	void
	show(
		std::ostream & out
	) const;

	friend
	std::ostream &
	operator<< (
		std::ostream & out,
		const Sqlite3ConnectionManager & sqlite3_connection_manager
	);



private:  // Private Functions

	void
	handle_db_errors(
		sqlite3 * db_connection
	) const ;


private:  // Private Data

	//TODO track open connections.
	//static std::map< char * , int > open_connection_counts_;
	ConnectionLists free_connections_;



	/// @brief singleton instance
	static Sqlite3ConnectionManagerOP instance_;


}; // Sqlite3ConnectionManager


} // utility
} // sqlite_database


#endif // DB_SQLITE3

#endif // INCLUDE_utility_sql_database_sqlite3_connection_manager_HH
