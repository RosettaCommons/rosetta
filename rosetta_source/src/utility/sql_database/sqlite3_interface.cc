// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   utility/sql_database/Sqlite3Interface.cc
/// @brief  Easy Access to sqlite3 database
/// @author Matthew O'Meara

//Unit Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/sqlite3_interface.hh>

#ifndef DB_SQLITE3

namespace utility {
namespace sql_database {

Sqlite3Interface::Sqlite3Interface() {}

Sqlite3Interface::Sqlite3Interface( Sqlite3Interface const & ) :
	utility::pointer::ReferenceCount()
{}

Sqlite3Interface::~Sqlite3Interface() {}

} // sql_database namespace
} // utility namespace

#else
// This requires the external dependency of the sqlite3 library
// To use compile with $scons.py extras=sqlite

// Package Header
#include <utility/sql_database/sqlite3_connection_manager.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <sqlite3.h>

//C++ Headers
#include <string>
#include <iostream>
#include <sstream>
namespace utility{
namespace sql_database{

using namespace std;
using namespace core;
using namespace utility;

	Sqlite3Interface::Sqlite3Interface() :
		sqlite3_connection_manager_( Sqlite3ConnectionManager::get_instance() ),
		database_fname_(),
		initialized_(false),
		db_( NULL ),
		active_column_( 0 ),
		active_column_total_( 0 ),
		active_stmt_( NULL ),
		column_count_()
	{}

	Sqlite3Interface::Sqlite3Interface(
		string const & database_fname ) :
		sqlite3_connection_manager_( Sqlite3ConnectionManager::get_instance() ),
		database_fname_( database_fname ),
		initialized_(false),
		db_( NULL ),
		active_column_( 0 ),
		active_column_total_( 0 ),
		active_stmt_( NULL ),
		column_count_()
	{
		open_connection(database_fname_);
	}

	Sqlite3Interface::Sqlite3Interface(
		Sqlite3Interface const & src ) :
		ReferenceCount( src ),
		sqlite3_connection_manager_( src.sqlite3_connection_manager_ ),
		database_fname_( src.database_fname_ ),
		initialized_(),
		db_( src.db_ ),
		active_column_( src.active_column_ ),
		active_column_total_( src.active_column_total_ ),
		active_stmt_( src.active_stmt_ ),
		column_count_(src.column_count_)
	{
		if (src.initialized_){
			open_connection(database_fname_);
			initialized_ = true;
		}
	}

	Sqlite3Interface::~Sqlite3Interface()
	{
		close_connection();
	}

	void
	Sqlite3Interface::open_connection(
		string const & database_fname
	) {
		if(initialized_){
			utility_exit_with_message("Attempting to a open connection to database '" +
																database_fname + "'while a connection to database '" +
																database_fname_ + "' is already open." );
		}
		database_fname_ = database_fname;
		sqlite3_connection_manager_->get_connection(database_fname_.c_str(), db_);
		initialized_ = true;
	}

	void
	Sqlite3Interface::close_connection(
	){
		if(initialized_){
			sqlite3_connection_manager_->free_connection(database_fname_.c_str(), db_);
			initialized_ = false;
		}

	}

	Sqlite3Interface::Sqlite3Const
	Sqlite3Interface::begin_transaction() {
		string const sql_stmt( "BEGIN DEFERRED TRANSACTION;" );
		Sqlite3Const ret = execute_sql(sql_stmt);
		return ret;
	}

	Sqlite3Interface::Sqlite3Const
	Sqlite3Interface::end_transaction() {
		string const sql_stmt( "END TRANSACTION;" );
		Sqlite3Const ret = execute_sql(sql_stmt);
		return ret;
	}

	string
	Sqlite3Interface::filename() const {
		return database_fname_;
	}

//  // Do not use by default for backwards compatibility
//	Size
//	Sqlite3Interface::memory_used(){
//		sqlite3_memory_highwater(1);
//		return 0;
//	}

	Sqlite3Interface::Sqlite3Const
	Sqlite3Interface::execute_sql(
		string const & sql_stmt
	){
		assert( initialized_ );
		Sqlite3Const ret = sqlite3_exec(db_, sql_stmt.c_str(),0,0,0);
		handle_db_errors( sql_stmt );

		return ret;
	}

	bool
	Sqlite3Interface::table_exists(
		string const & table_name
	){
		assert( initialized_ );
		string sql_stmt = "SELECT name FROM sqlite_master WHERE type='table' AND name='"+table_name+"';";
		sqlite3_stmt* stmt;
		sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &stmt, 0);
		handle_db_errors(stmt);
		bool exists(sqlite3_step(stmt) == SQLITE_ROW);
		sqlite3_finalize(stmt);
		handle_db_errors(stmt);
		return exists;
	}

	Size
	Sqlite3Interface::get_column_count(
		string const & tbl_name
	){

		assert( initialized_ );

		// This 'pragma' returns one row per column in table.
		// See: http://www.sqlite.org/pragma.html#pragma_table_info (accessed 30/7/10)
		string sql_stmt = "pragma table_info("+ tbl_name + ")";

		sqlite3_stmt * stmt;
		Sqlite3Const ret =
			sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &stmt, 0);
		handle_db_errors(stmt);
		ret = sqlite3_step(stmt);
		int n_columns = 0;
		while (ret == SQLITE_ROW){
			n_columns++;
			ret = sqlite3_step(stmt);
		}
		sqlite3_finalize(stmt);
		handle_db_errors(stmt);
		return n_columns;
	}

	Sqlite3Interface &
	Sqlite3Interface::operator<<(const begin_row & br ){
		assert( initialized_ );

		if( column_count_.find(br.t_name_) == column_count_.end() ){
			if( !table_exists( br.t_name_ ) ){
				cout << "sqlit3_interface ERROR: Attempting to insert into a table '" << br.t_name_ << "' but it doesn't exist in the database.  To create a table execute a 'CREATE TABLE statement on the database." << std::endl;
				utility_exit();
			}
			column_count_[br.t_name_] = get_column_count( br.t_name_ );
		}

		active_column_ = 1;
		active_column_total_ = column_count_[br.t_name_];

		if( active_column_total_ == 0){
			cout << "sqlite3_interface ERROR: Attempting to insert row in to table '" << br.t_name_ << "', but it has no columns." << std::endl;
			utility_exit();
		}

		stringstream prepare_stmt;
		prepare_stmt << "INSERT into " << br.t_name_ << " VALUES (";
		for(Size i=1; i<=active_column_total_; i++){
			prepare_stmt << "?" << i << (i< active_column_total_ ? "," : ");" );
		}

		sqlite3_prepare_v2(db_, prepare_stmt.str().c_str(), -1, &active_stmt_, 0);
		handle_db_errors( active_stmt_ );
		return *this;
	}

	// applicator of a parameterless manipulator
	// see for example:
	// http://appresario.com/online_books/programming_books/c++_practical_programming/c++_practical_programming_093.html
	Sqlite3Interface &
	Sqlite3Interface::operator<< (Sqlite3Interface &(*pf)(Sqlite3Interface &)){
		return pf(*this);
	}

	// bind null value to prepare statement
	//
	Sqlite3Interface &
	Sqlite3Interface::sqlite3_null(Sqlite3Interface& si){

		assert( si.initialized_ );

		if( si.active_stmt_ ){
			sqlite3_bind_null(si.active_stmt_, si.active_column_);
			si.handle_db_errors( si.active_stmt_ );
			si.active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}
		return si;
	}


	Sqlite3Interface &
	Sqlite3Interface::end_row(Sqlite3Interface & si){

		assert( si.initialized_ );

		if( si.active_stmt_ ){

			assert( si.active_column_ -1 == si.active_column_total_ );

			Sqlite3Interface::Sqlite3Const ret = sqlite3_step(si.active_stmt_);
			if( ret != SQLITE_DONE ){
				si.handle_db_errors( si.active_stmt_ );
			}
			sqlite3_finalize(si.active_stmt_);
			si.handle_db_errors( si.active_stmt_ );
		} else {
			cout << "Attempting to end a row when ther is no row statement to finish!" << endl;
		}
		return si;
	}

	Sqlite3Interface &
	Sqlite3Interface::operator<< ( int const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_int(active_stmt_, active_column_, value);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}

	Sqlite3Interface &
	Sqlite3Interface::operator<< ( Size const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_int(active_stmt_, active_column_, value);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}


	Sqlite3Interface &
	Sqlite3Interface::operator<< ( Real const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_double(active_stmt_, active_column_, value);
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}

	Sqlite3Interface &
	Sqlite3Interface::operator<< ( string const & value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_text(active_stmt_, active_column_, value.c_str(), -1, SQLITE_TRANSIENT);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}

	Sqlite3Interface &
	Sqlite3Interface::operator<<( char const * const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_text(active_stmt_, active_column_, value, -1, SQLITE_TRANSIENT);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}



	Sqlite3Interface &
	Sqlite3Interface::operator<<( char const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_text(active_stmt_, active_column_, &value, 1, SQLITE_TRANSIENT);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}



	Sqlite3Interface &
	Sqlite3Interface::operator<< ( bool const value ){

		assert( initialized_ );

		if( active_stmt_ ){
			sqlite3_bind_int(active_stmt_, active_column_, value);
			handle_db_errors( active_stmt_ );
			active_column_++;
		} else {
			cout << "Attempting to push a value but no statement to bind it to!" << endl;
		}

		return *this;
	}

	Size
	Sqlite3Interface::last_key() const {

		assert( initialized_ );

		Size key = sqlite3_last_insert_rowid(db_);
		handle_db_errors( active_stmt_ );
		return key;
	}

	void
	Sqlite3Interface::handle_db_errors(
		sqlite3_stmt * stmt
	) const {

		assert( initialized_ );

		if( sqlite3_errcode(db_) ){
			string errstr = sqlite3_errmsg(db_);
			cout << "SQLITE ERROR MESSAGE: " << errstr << endl;
			if ( stmt ){
				cout << "While executing: '" << sqlite3_sql(active_stmt_) << "'" << endl;
			}
			utility_exit();
		}
	}

	void
	Sqlite3Interface::handle_db_errors(
		string const & sql_stmt
	) const {

		assert( initialized_ );

		if( sqlite3_errcode(db_) ){
			string errstr = sqlite3_errmsg(db_);
			cout << "SQLITE ERROR MESSAGE: " << errstr << endl;
			if ( sql_stmt.compare("") != 0 ){
				cout << "While executing: '" << sql_stmt << "'" << endl;
			}
			utility_exit();
		}
	}



	void
	Sqlite3Interface::show(
		ostream & out
	) const {
		out << "Sqlite3Interface (For Debugging):" << endl;
		out << "\tthis : '" << this << "'" << endl;
		out << "\tdatabase_fname_ : '" << database_fname_ << "'" << endl;
		out << "\tinitialized_ : '" << initialized_ << "'" << endl;
		out << "\tdb_ : '" << db_ << "'" << endl;
		out << "\tactive_column_ : '" << active_column_ << "'" << endl;
		out << "\tactive_column_total_ : '" << active_column_total_ << "'" << endl;
		out << "\tactive_stmt_ : '" << active_stmt_ << "'" << endl;
	}

	ostream &
	operator<< (
		ostream & out,
		const Sqlite3Interface & sqlite3_interface
	) {
		sqlite3_interface.show( out );
		return out;
	}




} // sql_database
} // utility

#endif // DB_SQLITE3
