// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/sqlite3_connection_manager.cc
/// @brief  Easy Access to sqlite3 database; manage connections
/// @author Matthew O'Meara

#ifdef DB_SQLITE3

#include <utility/sql_database/sqlite3_connection_manager.hh>
#include <utility/exit.hh>

#include <sqlite3.h>

#include <string>
#include <iostream>
#include <list>

using std::string;
using std::endl;
using std::cout;
using std::list;
using std::ostream;


namespace utility{
namespace sql_database{

	//Initialize static member variable
	Sqlite3ConnectionManagerOP Sqlite3ConnectionManager::instance_;

	Sqlite3ConnectionManager::ConnectionList::ConnectionList() :
		connection_list_()
	{}

	Sqlite3ConnectionManager::ConnectionList::ConnectionList(
		sqlite3 * & db_connection
	) :
		connection_list_( 1, db_connection )
	{}

	Sqlite3ConnectionManager::ConnectionList::~ConnectionList(
	) {}

	void
	Sqlite3ConnectionManager::ConnectionList::free_all(){
		for( list< sqlite3 * >::iterator
			c_it = connection_list_.begin(),
			c_end = connection_list_.end();
			c_it != c_end; ++c_it ){
			sqlite3_close( *c_it );
		}
	}

	void
	Sqlite3ConnectionManager::ConnectionList::push_front(
		sqlite3 * & db_connection
	) {
		connection_list_.push_front( db_connection );
	}

	void
	Sqlite3ConnectionManager::ConnectionList::extract_front(
		sqlite3 * & db_connection
	) {
		if( connection_list_.empty() ){
			db_connection = NULL;
		} else {
			db_connection = connection_list_.front();
			connection_list_.pop_front();
		}
	}

	bool
	Sqlite3ConnectionManager::ConnectionList::empty(
	) {
		return connection_list_.empty();
	}

	void
	Sqlite3ConnectionManager::ConnectionList::show(
		ostream & out
	) const {
		for( std::list< sqlite3 * >::const_iterator
			db_con_it = connection_list_.begin(),
			db_con_end = connection_list_.end();
			db_con_it != db_con_end; ++db_con_it ) {
			out << "\tsqlite3 * : '" << *db_con_it << "'" << endl;
		}
	}

	ostream &
	operator<< (
		ostream & out,
		const Sqlite3ConnectionManager::ConnectionList & connection_list
	) {
		connection_list.show( out );
		return out;
	}


	Sqlite3ConnectionManagerOP
	Sqlite3ConnectionManager::get_instance() {
		if ( !instance_.get() ) {
			instance_ = new Sqlite3ConnectionManager();
		}

		return instance_;
	}

	Sqlite3ConnectionManager::Sqlite3ConnectionManager() {}


	Sqlite3ConnectionManager::~Sqlite3ConnectionManager() {
		// Assume all open connections have been freed.
		// All free connections are closed on destruction.
	}

	int
	Sqlite3ConnectionManager::get_connection(
		string const & db_fname,
		sqlite3 * & db_connection
	){
		int ret = 0; // SQLITE_OK

		ConnectionListOP & cl( free_connections_[db_fname] );
		if( !cl || cl->empty() ){
			cl = new ConnectionList();
			int ret = sqlite3_open(db_fname.c_str(), &db_connection);
			if(!ret){
				cout << "Opened connection to sqlite database '" << db_fname << "'" << endl;
			}
		} else {
			cl->extract_front( db_connection );
		}
		return ret;
	}

	void
	Sqlite3ConnectionManager::free_connection(
		string const db_fname,
		sqlite3 * & db_connection
	){
		ConnectionListOP cl( free_connections_[ db_fname ] );
		if (!cl){
			cl = new ConnectionList( db_connection );
		} else {
			cl->push_front( db_connection );
		}
	}

	void
	Sqlite3ConnectionManager::show(
		ostream & out
	) const {
		cout << "\tthis : '" << this << "'" << endl;
		for( ConnectionLists::const_iterator
			cl_it = free_connections_.begin(),
			cl_end = free_connections_.end();
			cl_it != cl_end; ++cl_it ){
			string key_( cl_it->first );
			out << "\tDatabase : '" << key_ << "'" << endl;
			out << "\t" << *(cl_it->second);
		}
	}

	ostream &
	operator<< (
		ostream & out,
		const Sqlite3ConnectionManager & sqlite3_connection_manager
	) {
		sqlite3_connection_manager.show( out );
		return out;
	}


	void
	Sqlite3ConnectionManager::handle_db_errors(
		sqlite3 * db_connection
	) const {
		if( sqlite3_errcode(db_connection) ){
			string errstr = sqlite3_errmsg(db_connection);
			cout << "SQLITE ERROR MESSAGE: " << errstr << endl;
			utility_exit();
		}
	}

} // sql_database
} // utility


#endif // DB_SQLITE3
