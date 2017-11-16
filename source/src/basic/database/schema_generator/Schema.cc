// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/schema_generator/Schema.cc
///
/// @brief Construct a database backend independant schema
/// @author Tim Jacobs

//Unit
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Index.hh>

#include <basic/database/sql_utils.hh>

#include <basic/mpi/MessageListenerFactory.hh>
#include <basic/mpi/util.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <platform/types.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/mpi_util.hh>
#include <utility/sql_database/types.hh>

#include <string>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <cctype>

// External
#include <cppdb/frontend.h>

namespace basic {
namespace database {
namespace schema_generator {

using platform::Size;
using std::string;
using std::endl;
using std::stringstream;
using utility::sql_database::sessionOP;
using basic::mpi::MessageListenerOP;
using basic::mpi::MessageListenerFactory;
using basic::mpi::request_data_from_head_node;
using basic::mpi::send_data_to_head_node;
using basic::database::table_exists;
using utility::mpi_rank;
using cppdb::statement;

static basic::Tracer TR( "basic.database.schema_generator.Schema" );


Schema::Schema(std::string table_name):
	table_name_(table_name)
{
	init();
}

Schema::Schema(std::string table_name, PrimaryKey primary_key):
	table_name_(table_name),
	primary_key_(primary_key)
{
	init();
}

Schema::Schema(
	Schema const & src
) :
	primary_key_(src.primary_key_),
	columns_(src.columns_),
	foreign_keys_(src.foreign_keys_),
	constraints_(src.constraints_),
	indices_(src.indices_)
{}

void
Schema::init(){
	//Add primary key columns to schema list
	Columns key_columns = primary_key_.columns();
	this->columns_.insert( columns_.end(), key_columns.begin(), key_columns.end() );

	// Table names should all be lower case
	std::transform(
		table_name_.begin(), table_name_.end(), table_name_.begin(),
		(int(*)(int)) std::tolower);
}

void
Schema::add_foreign_key(
	ForeignKey key
){
	this->foreign_keys_.push_back(key);
	//if the foreign key is also a primary key it will have already been added

	Columns key_cols = key.columns();

	for ( Size i=1; i <= key_cols.size(); ++i ) {
		if ( !this->columns_.contains(key_cols[i]) ) {
			this->columns_.push_back(key_cols[i]);
		}
	}
}

void
Schema::add_column(
	Column column
){
	//Don't add a column more than once
	if ( !this->columns_.contains(column) ) {
		this->columns_.push_back(column);
	}
}

void
Schema::add_constraint(
	ConstraintOP constraint
){
	this->constraints_.push_back(constraint);
}

void
Schema::add_index(
	Index index
) {
	indices_.push_back(index);
}

std::string Schema::print( sessionOP db_session ) const
{
	return table_schema_statements(db_session) + table_init_statements(db_session);
}

std::string Schema::table_schema_statements( sessionOP db_session ) const
{
	stringstream schema_string;
	schema_string << "CREATE TABLE IF NOT EXISTS " << table_name_ << "(\n\t";

	for ( Columns::const_iterator it=columns_.begin(), end = columns_.end(); it != end; ++it ) {
		if ( it!=columns_.begin() ) {
			schema_string << ",\n\t";
		}
		schema_string << it->print(db_session);
	}

	for ( size_t i=1; i<=foreign_keys_.size(); i++ ) {
		schema_string << ",\n\t" << foreign_keys_[i].print(db_session);
	}

	Columns const & keys(primary_key_.columns());

	if ( keys.size() > 0 ) {
		switch(db_session->get_db_mode()) {
		case utility::sql_database::DatabaseMode::mysql:
		case utility::sql_database::DatabaseMode::postgres :
			schema_string << ",\n\t" << primary_key_.print(db_session);
			break;
		case utility::sql_database::DatabaseMode::sqlite3 :
			//Prevent adding the primary key twice - this will happen if you have an autoincrementing primary key in sqlite3

			if ( !(keys.size()==1 && keys.begin()->auto_increment()) ) {
				schema_string << ",\n\t" << primary_key_.print(db_session);
			}
			break;
		default :
			utility_exit_with_message(
				"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
		}
	}

	for ( size_t i=1; i<=constraints_.size(); i++ ) {
		schema_string << ",\n\t" << constraints_[i]->print(db_session);
	}

	schema_string << ");\n";

	for ( Size i=1; i <= indices_.size(); ++i ) {
		schema_string << indices_[i].print(table_name_, db_session);
		schema_string << "\n";
	}

	return schema_string.str();
}

std::string Schema::table_init_statements( sessionOP db_session ) const
{
	stringstream init_string;
	for ( Columns::const_iterator it=columns_.begin(), end = columns_.end(); it != end; ++it ) {
		if ( it->auto_increment() && it->auto_increment_base() != 0 ) {
			switch(db_session->get_db_mode()){
			case utility::sql_database::DatabaseMode::postgres :
				// Postgresql creates an implicit sequence of the form <table_name>_<column_name>_seq to handle auto_increment values.
				//
				// http://www.postgresql.org/docs/9.2/static/datatype-numeric.html
				// see 8.1.4 Serial Types
				init_string << "ALTER SEQUENCE " << table_name_ << "_" << it->name() << "_seq" << " MINVALUE " << it->auto_increment_base() << ";\n";
				break;
			case utility::sql_database::DatabaseMode::mysql :
				// http://dev.mysql.com/doc/refman/5.6/en/example-auto-increment.html
				init_string << "ALTER TABLE " << table_name_ << " AUTO_INCREMENT = " << it->auto_increment_base() << ";\n";
				break;
			case utility::sql_database::DatabaseMode::sqlite3 :
				// Sqlite3 autoincrement is managed by the db table, sqlite_sequence. ROWID (and autoincrement value) are monotonically
				// increased from from current value in sqlite_sequence. In the trivial case, the next value is sqlite_sequence + 1.
				//
				// http://www.sqlite.org/autoinc.html
				init_string << "INSERT INTO sqlite_sequence SELECT \"" << table_name_ << "\", " << it->auto_increment_base() - 1 << "  WHERE NOT EXISTS (SELECT 1 FROM sqlite_sequence WHERE NAME = \"" << table_name_ << "\");\n";
				break;
			default :
				utility_exit_with_message(
					"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
			}
		}
	}

	return init_string.str();
}

//Write this schema to the database
void Schema::write(sessionOP db_session)
{
	std::string schema_statement = table_schema_statements(db_session);
	std::string init_statements = table_init_statements(db_session);

	try
{
		try
{
			/*
			* TODO alexford This abort logic still allows race conditions if writes are executed from multiple frontends
			* This should be resolved by adding a meta-table to track schema initialization protected by transactions, or
			* by running all table declations
			*/
			check_table_and_perform_write(db_session, schema_statement, init_statements);
		}
catch(cppdb::cppdb_error & except)
{
	TR.Debug << "Error writing schema, retrying:\n" << except.what() << std::endl;
	check_table_and_perform_write(db_session, schema_statement, init_statements);
}
	}
catch(cppdb::cppdb_error & except)
{
	TR.Error
		<< "ERROR writing schema after retry.\n"
		<< print(db_session) << std::endl;
	TR.Error << except.what() << std::endl;
	utility_exit();
}
}

void Schema::check_table_and_perform_write(
	utility::sql_database::sessionOP db_session,
	std::string schema_statement,
	std::string init_statements) const
{
	cppdb::transaction guard(*db_session);
	statement stmt;

	if ( init_statements != "" && table_exists(db_session, table_name_) ) {
		TR.Debug << "Table with init statement exists, skipping declaration: " << table_name_ << std::endl;
		return;
	}

	if ( !table_exists(db_session, table_name_) ) {
		// kylebarlow - We don't need to run a "CREATE TABLE IF NOT EXISTS" query if the table already exists
		//   Running lots of those queries results in problems waiting for table metadata locks
		//   This problem remains unfixed in MySQL as of version 5.5.32
		TR.Debug << "Writing schema for table: " << table_name_ << std::endl;
		TR.Trace << schema_statement << std::endl;

		stmt = (*db_session) << schema_statement;
		stmt.exec();
	}

	if ( init_statements != "" ) {
		TR.Debug << "Writing init for table: " << table_name_ << std::endl;
		TR.Trace << init_statements << std::endl;

		stmt = (*db_session) << init_statements;
		stmt.exec();
	}

	guard.commit();
}

} // schema_generator
} // namespace database
} // namespace utility
