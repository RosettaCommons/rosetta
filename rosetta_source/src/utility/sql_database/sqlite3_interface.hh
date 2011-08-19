// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/sqlite3_interface.hh
/// @brief  Easy Access to sqlite3 database
/// @author Matthew O'Meara



#ifndef INCLUDE_utility_sql_database_sqlite3_interface_HH
#define INCLUDE_utility_sql_database_sqlite3_interface_HH

// Unit Headers
#include <utility/sql_database/sqlite3_interface.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifndef DB_SQLITE3
// this interface requires the DB_SQLITE3 compilation flag
// Enable with $./scons.py extras=sqlite

namespace utility{
namespace sql_database{

// define the class even if DB_SQLITE3 is not used because it's type
// information may be needed for for classes that optionally take a
// database interface.  In that case, however, use the 'null_interface'
class Sqlite3Interface : public utility::pointer::ReferenceCount {
public:
	Sqlite3Interface();

	Sqlite3Interface( Sqlite3Interface const & );

	virtual
	~Sqlite3Interface();

};

static Sqlite3Interface null_interface;

} // sql_database namespace
} // utility namespace

#else // DB_SQLITE3

#include <utility/sql_database/sqlite3_connection_manager.fwd.hh>

// Core Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


// C++ Headers
#include <string>
#include <map>

// External Headers
#include <sqlite3.h>

namespace utility {
namespace sql_database {



class Sqlite3Interface : public utility::pointer::ReferenceCount {

public:

	class begin_row {

	public:
		std::string const t_name_;

		begin_row(std::string const tbl_name ) :
			t_name_( tbl_name )
		{}

	};

	typedef int Sqlite3Const;

	Sqlite3Interface();

	Sqlite3Interface(
		std::string const & database_fname
	);

	Sqlite3Interface( Sqlite3Interface const & src );

	virtual ~Sqlite3Interface();

	void
	open_connection(
		std::string const & database_fname
	);

	void
	close_connection();

	Sqlite3Const
	begin_transaction();

	Sqlite3Const
	end_transaction();

//  // Do not use by default for backwards compatibility
//	core::Size
//	memory_used();
	bool
	table_exists( std::string const & table_name );

	std::string
	filename() const;

	Sqlite3Const
	execute_sql( std::string const & sql_stmt );

//	void
//	get_table(
//						std::string const & selet_stmt,
//						vector1<

	core::Size
	get_column_count(
		std::string const & tbl_name);

	Sqlite3Interface &
	operator<<(const begin_row & br );

	Sqlite3Interface &
	operator<< (Sqlite3Interface &(*pf)(Sqlite3Interface &));

	// manipulator
	static
	Sqlite3Interface &
	sqlite3_null (Sqlite3Interface & si);

	// manipulator
	static
	Sqlite3Interface &
	end_row(Sqlite3Interface & si);

	Sqlite3Interface &
	operator <<( int const value );

	Sqlite3Interface &
	operator <<( core::Size const value );

	Sqlite3Interface &
	operator <<( core::Real const value );

	Sqlite3Interface &
	operator <<( std::string const & value );

	Sqlite3Interface &
	operator <<( char const * const value );

	Sqlite3Interface &
	operator <<( char const value );

	Sqlite3Interface &
	operator <<( bool const value );

	core::Size
	last_key() const;


	void
	show(
		std::ostream & out
	) const;

	friend
	std::ostream &
	operator<< (
		std::ostream & out,
		const Sqlite3Interface & sqlite3_interface
	);



private:

	void
	handle_db_errors(
		sqlite3_stmt * stmt
	) const;

	void
	handle_db_errors(
		std::string const & stmt
	) const;


	Sqlite3ConnectionManagerOP sqlite3_connection_manager_;
	std::string database_fname_;
	bool initialized_;

	sqlite3* db_;

	core::Size active_column_;
	core::Size active_column_total_;
	sqlite3_stmt* active_stmt_;


	std::map< std::string const, core::Size > column_count_;

};

} // utilitiy
} // sql_database

#endif // DB_SQLITE3

#endif // INCLUDE_utility_sql_database_sqlite3_interface_HH
