// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/database/sql_utils.cc
/// @brief Database utility functions
/// @author Sam DeLuca
/// @author Matthew O'Meara

#ifdef USEMPI
#include <mpi.h>
#endif

#include <basic/database/sql_utils.hh>

#include <basic/Tracer.hh>                                     // for Tracer
#include <basic/options/keys/inout.OptionKeys.gen.hh>          // for databa...
#include <basic/options/keys/out.OptionKeys.gen.hh>            // for all, db
#include <basic/options/option.hh>                             // for Option...
#include <basic/resource_manager/ResourceManager.hh>           // for Resour...
#include <basic/resource_manager/util.hh>                      // for get_re...

#include <cppdb/errors.h>                                      // for cppdb_...
#include <cppdb/frontend.h>                                    // for statement

#include <boost/algorithm/string/predicate.hpp>                // for istart...
#include <boost/foreach.hpp>                                   // for auto_a...
#include <boost/iterator/iterator_facade.hpp>                  // for operat...
#include <boost/lexical_cast.hpp>                              // for bad_le...
#include <boost/mpl/bool.hpp>                                  // for bool_
#include <boost/mpl/bool_fwd.hpp>                              // for false_
#include <boost/token_functions.hpp>                           // for char_s...
#include <boost/tokenizer.hpp>                                 // for tokenizer
#include <numeric/random/random.hh>                            // for rg
#include <platform/types.hh>                                   // for Size
#include <utility/excn/Exceptions.hh>                          // for EXCN_M...
#include <utility/exit.hh>                                     // for utilit...
#include <utility/file/PathName.hh>                            // for PathName
#include <utility/options/BooleanOption.hh>                    // for Boolea...
#include <utility/options/IntegerOption.hh>                    // for Intege...
#include <utility/options/PathOption.hh>                       // for PathOp...
#include <utility/options/ScalarOption_T_.hh>                  // for Scalar...
#include <utility/options/StringOption.hh>                     // for String...
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>  // for sessionOP
#include <utility/sql_database/DatabaseSessionManager.hh>      // for session
#include <utility/string_util.hh>                              // for join
#include <utility/vector1.hh>                                  // for vector1
#include <utility/tag/Tag.hh>                                  // for Tag

#include <cstddef>                                             // for size_t
#include <iosfwd>                                              // for string
#include <memory>                                              // for shared...
#include <sstream>                                             // for operat...
#include <string>                                              // for basic_...
#include <vector>                                              // for vector

#ifndef WIN32
#include <unistd.h>                                            // for usleep
#endif


#ifdef WIN32
#include <windows.h>
#endif

using std::string;
using std::stringstream;
using utility::sql_database::sessionOP;
using platform::Size;
using cppdb::statement;
using cppdb::cppdb_error;
using cppdb::result;
using utility::vector1;
using namespace utility::sql_database;

namespace basic {
namespace database {

static THREAD_LOCAL basic::Tracer TR( "basic.database.sql_utils" );


sessionOP
get_db_session() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		option[inout::dbms::database_name],
		utility::sql_database::TransactionMode::standard,
		0,
		option[inout::dbms::pq_schema]);
}


sessionOP
get_db_session(
	string const & db_name,
	string const & pq_schema
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		database_mode_from_name(option[inout::dbms::mode]),
		utility::sql_database::TransactionMode::standard,
		0,
		db_name,
		pq_schema);
}

sessionOP
get_db_session(
	string const & db_name,
	TransactionMode::e transaction_mode,
	Size chunk_size,
	string const & pq_schema
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		database_mode_from_name(option[inout::dbms::mode]),
		transaction_mode,
		chunk_size,
		db_name,
		pq_schema);
}

utility::sql_database::sessionOP
get_db_session(
	utility::sql_database::DatabaseMode::e db_mode,
	std::string const & db_name,
	std::string const & pq_schema){

	return get_db_session(
		db_mode,
		utility::sql_database::TransactionMode::standard,
		0,
		db_name,
		pq_schema);
}

sessionOP
get_db_session(
	DatabaseMode::e db_mode,
	TransactionMode::e transaction_mode,
	Size chunk_size,
	string const & db_name,
	string const & pq_schema
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	switch(db_mode) {
	case DatabaseMode::sqlite3 :

		if (
				option[inout::dbms::host].user() ||
				option[inout::dbms::user].user() ||
				option[inout::dbms::password].user() ||
				option[inout::dbms::port].user() ) {
			utility_exit_with_message(
				"You have specified options for a client-server database "
				"but the database mode is sqlite3. "
				"Please specify -inout:dbms:mode <db_mode>.");
		}

		if ( pq_schema.compare("") ) {
			TR.Warning
				<< "You have specified a postgres schema but using a sqlite3 database. "
				<< "To use postgres, please specify -inout:dbms:mode postgres"
				<< std::endl;
		}

		return DatabaseSessionManager::get_instance()->get_session_sqlite3(
			db_name,
			transaction_mode,
			chunk_size,
			option[inout::dbms::readonly],
			db_partition_from_options(db_mode));

	case DatabaseMode::mysql :

		if ( option[inout::dbms::readonly] ) {
			utility_exit_with_message(
				"Restricting access to a mysql database is done at the user level "
				"rather that the connection level. "
				"So requesting a readonly connection cannot fullfilled.");
		}

		if ( pq_schema.compare("") ) {
			TR.Warning
				<< "You have specified a postgres schema but using a mysql database. "
				<< "To use postgres, please specify -inout:dbms:mode postgres"
				<< std::endl;
		}

		if ( !(
				option[inout::dbms::host].user() &&
				option[inout::dbms::user].user() &&
				option[inout::dbms::password].user() &&
				option[inout::dbms::port].user()) ) {
			utility_exit_with_message(
				"To connect to a mysql database you must specify "
				"-inout:dbms:host -inout:dbms:user -inout:dbms:password and "
				"-inout:dbms:port");
		}

		// Call to get partition performs option validation.
		db_partition_from_options(db_mode);

		return DatabaseSessionManager::get_instance()->get_session_mysql(
			db_name,
			transaction_mode,
			chunk_size,
			option[inout::dbms::host],
			option[inout::dbms::user],
			option[inout::dbms::password],
			option[inout::dbms::port]);


	case DatabaseMode::postgres :

		if ( option[inout::dbms::readonly] ) {
			utility_exit_with_message(
				"Restricting access to a postgres database is done at the user level "
				"rather that the connection level. So requesting a readonly connection "
				"cannot fullfilled.");
		}


		if ( !(
				option[inout::dbms::host].user() &&
				option[inout::dbms::user].user() &&
				option[inout::dbms::password].user() &&
				option[inout::dbms::port].user()) ) {
			utility_exit_with_message(
				"To connect to a postgres database you must specify "
				"-inout:dbms:host -inout:dbms:user -inout:dbms:password "
				"and -inout:dbms:port");
		}

		// Call to get partition performs option validation.
		db_partition_from_options(db_mode);

		return DatabaseSessionManager::get_instance()->get_session_postgres(
			db_name,
			transaction_mode,
			chunk_size,
			pq_schema,
			option[inout::dbms::host],
			option[inout::dbms::user],
			option[inout::dbms::password],
			option[inout::dbms::port]);
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_mode) + "'");
	}
	return 0;
}

/// @
platform::SSize db_partition_from_options( DatabaseMode::e db_mode )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[inout::dbms::separate_db_per_mpi_process] && option[inout::dbms::database_partition].user() ) {
		utility_exit_with_message(
			"The -inout:dbms:separate_db_per_mpi_process and -inout::dbms::database_partition options are mutually exclusive.");
	} else if ( option[inout::dbms::database_partition].user() ) {
		if ( option[inout::dbms::database_partition] < 0 ) {
			stringstream error;
			error << "Invalid -inout::dbms::database_partition specified, must be positive integer value: " << option[inout::dbms::database_partition] << std::endl;
			utility_exit_with_message(error.str());
		}
	}

	switch(db_mode)
			{
			case DatabaseMode::sqlite3 :

				if ( !option[inout::dbms::database_partition].user() ) {
					return resolve_db_partition(
						option[inout::dbms::separate_db_per_mpi_process]);
				} else {
					return resolve_db_partition(
						option[inout::dbms::separate_db_per_mpi_process],
						option[inout::dbms::database_partition]);
				}

				break;

			case DatabaseMode::mysql:
			case DatabaseMode::postgres :
				if ( option[inout::dbms::separate_db_per_mpi_process] || option[inout::dbms::database_partition].user() ) {
					utility_exit_with_message(
						"The -inout:dbms:separate_db_per_mpi_process and -inout::dbms::database_partition flags only apply to "
						"sqlite3 databases.");
				}

				break;
			}

	return -1;
}

platform::SSize resolve_db_partition(
	bool partition_by_mpi_process,
	platform::SSize manual_partition)
{
	if ( partition_by_mpi_process ) {
#ifdef USEMPI
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		return mpi_rank;
#endif
	}

	if ( manual_partition >= 0 ) {
		return manual_partition;
	}

	return -1;
}

statement safely_prepare_statement(
	string const & statement_string,
	sessionOP db_session)
{
	statement stmt;
	try
{
		stmt = db_session->prepare(statement_string);
		return stmt;
	}catch(cppdb_error const & error)
{
		TR.Error << " Failed to safely prepare the following statement: " << std::endl;
		TR.Error << statement_string << std::endl;
		TR.Error << error.what() << std::endl;

		utility_exit_with_message(error.what());
	}
	return stmt; //there's no way this should happen
}

void
safely_write_to_database(
	statement & statement
) {

	platform::Size cycle = 0;
	while ( true )
			{
		try
{
			statement.exec();
			return;
		}catch(cppdb::bad_value_cast & except)
{

			utility_exit_with_message(except.what());
		}catch(cppdb::empty_row_access & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_column & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_placeholder & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::multiple_rows_query & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::not_supported_by_backend & except)
{
			utility_exit_with_message(except.what());
		} catch(cppdb::null_value_fetch & except)
{
			utility_exit_with_message(except.what());
		} catch(cppdb::cppdb_error & except)
{
#ifdef USEMPI
			if( std::string(except.what()) == "database is locked"){
				stringstream err_msg;
				err_msg
					<< "database is locked" << std::endl
					<< std::endl
					<< "PSA: If this is an sqlite3 session, running under MPI, and you are using the database primarily for writing output," << std::endl
					<< "consider using the 'separate_db_per_mpi_process' option if you're not already using it." << std::endl
					<< "To do this add separate_db_per_mpi_process=1 to a RosettaScripts/resource_definitions xml tag that takes database options" << std::endl
					<< "or add -inout:dbms:separate_db_per_mpi_process to the command line or flags file" << std::endl
					<< "This option will append '_<mpi_rank>' to the database filename for each mpi process. Once the run has finished," << std::endl
					<< "these databases can be merged together using " << std::endl
					<< std::endl
					<< "        bash /path/to/rosetta_tests/features/sample_sources/merge.sh <output_db> <input_db_part_1> [<input_db_part_2> [ ... ] ]" << std::endl
					<< std::endl
					<< "For more information see the database connection options page in the Rosetta Wiki: RosettaScripts_database_connection_options" << std::endl
					<< except.what() << std::endl;
				utility_exit_with_message(err_msg.str());
			}
#endif

			utility_exit_with_message(except.what());
		}

		cycle++;
	}
}

result
safely_read_from_database(
	statement & statement
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	platform::Size retry_limit = 10;
	platform::Size cycle = 1;
	while ( true )
			{
		try
{
			return statement.query();
		}catch(cppdb::bad_value_cast & except)
{

			utility_exit_with_message(except.what());
		}catch(cppdb::empty_row_access & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_column & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_placeholder & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::multiple_rows_query & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::not_supported_by_backend & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::null_value_fetch & except)
{
			utility_exit_with_message(except.what());
		}catch(cppdb::cppdb_error & except)
{
			if ( option[inout::dbms::retry_failed_reads] && cycle <= retry_limit ) {
				TR << "Caught an exception on db read: "  << except.what() << "\n" <<
					"pausing before another attempted read. This is attempt " << cycle << " of " << retry_limit <<std::endl;
#ifdef WIN32
				Sleep(1000);
#else
				//Sleep some amount between 100-2000 ms
				usleep(100+1900*numeric::random::rg().uniform());
#endif
			} else {
				utility_exit_with_message(except.what());
			}
		}
		cycle++;
	}
}


bool
table_exists(
	sessionOP db_session,
	string const & table_name
) {

	// TODO: handle when the current database is not the one from the
	// option system, can someone with mysql try this and see if it
	// works?
	// "SHOW TABLES IN database();"
	string statement_string;
	statement stmt;
	Size i(1);
	switch(db_session->get_db_mode()){
	case DatabaseMode::sqlite3 :
		statement_string = "SELECT name FROM sqlite_master WHERE name=?;";
		stmt = safely_prepare_statement(statement_string,db_session);
		break;
	case DatabaseMode::mysql :
		statement_string = "SHOW TABLES WHERE Tables_in_"+db_session->get_db_name()+" = ?;";
		stmt = safely_prepare_statement(statement_string,db_session);
		break;
	case DatabaseMode::postgres :
		if ( db_session->get_pq_schema() == "" ) {
			statement_string =
				"SELECT tablename \n"
				"FROM pg_catalog.pg_tables \n"
				"WHERE tablename = ?;";
			stmt = safely_prepare_statement(statement_string, db_session);
		} else {
			statement_string =
				"SELECT tablename FROM pg_catalog.pg_tables \n"
				"WHERE schemaname = ? AND tablename = ?;";
			stmt = safely_prepare_statement(statement_string, db_session);
			stmt.bind(i, db_session->get_pq_schema());
			i++;
		}
		break;
	default :
		utility_exit_with_message("unknown database mode");
	}

	stmt.bind(i, table_name);
	result res(stmt.query());

	return res.next();
}

//Simply (probably overly so) protection from SQL injection
void
check_statement_sanity(
	string sql
) {
	if ( !boost::istarts_with(sql, "SELECT") ) {
		utility_exit_with_message("ERROR: Database select statement is not safe! Only SELECT statements are allowed.");
	}

	int semicolon_count=0;
	for ( size_t i = 0; i < sql.size(); i++ ) {
		if ( sql[i] == ';' ) {
			semicolon_count++;
		}
	}
	if ( semicolon_count > 1 ) {
		utility_exit_with_message("Database select statement is not safe! Only 1 SQL statement is allowed");
	}
}

//This should ideally only be used for reference tables that have static data that needs to only be written once(ex: dssp_codes)
void
insert_or_ignore(
	string table_name,
	std::vector<string> column_names,
	std::vector<string> values,
	sessionOP db_session
){

	string statement_string="";
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::mysql : {
		statement_string = "INSERT IGNORE into " + table_name + " (";
		for ( size_t i=0; i<column_names.size(); i++ ) {
			statement_string += column_names[i];
			if ( i != column_names.size()-1 ) {
				statement_string += ",";
			}
		}

		statement_string+=") VALUES(";
		for ( size_t i=0; i<values.size(); i++ ) {
			statement_string += "?";
			if ( i != column_names.size()-1 ) {
				statement_string += ",";
			}
		}

		statement_string += ");";

		statement stmt( safely_prepare_statement(statement_string, db_session) );
		for ( size_t i=0; i<values.size(); i++ ) {
			stmt.bind(i+1, values[i]);
		}

		safely_write_to_database(stmt);
		break;
	}
	case utility::sql_database::DatabaseMode::postgres : {
		//This is a dirty postgres hack and seems to be the easiest workaround for lack of INSERT IGNORE support in postgres
		string select_statement_string = "SELECT * FROM " + table_name + " WHERE ";
		for ( size_t i=0; i<column_names.size(); i++ ) {
			select_statement_string += column_names[i] + "=?";
			if ( i != column_names.size()-1 ) {
				select_statement_string += " AND ";
			}
		}
		select_statement_string+=";";

		statement select_stmt ( safely_prepare_statement(
			select_statement_string, db_session
			));

		for ( size_t i=0; i<values.size(); i++ ) {
			select_stmt.bind(i+1, values[i]);
		}

		result res = safely_read_from_database(select_stmt);

		if ( !res.next() ) {
			statement_string += "INSERT into " + table_name + "(";
			for ( size_t i=0; i<column_names.size(); i++ ) {
				statement_string += column_names[i];
				if ( i != column_names.size()-1 ) {
					statement_string += ",";
				}
			}

			statement_string += ") VALUES(";
			for ( size_t i=0; i<values.size(); i++ ) {
				statement_string += "?";
				if ( i != column_names.size()-1 ) {
					statement_string += ",";
				}
			}

			statement_string += ");";

			statement stmt ( safely_prepare_statement(
				statement_string, db_session
				));

			for ( size_t i=0; i<values.size(); i++ ) {
				stmt.bind(i+1, values[i]);
			}

			safely_write_to_database(stmt);
		}
		break;
	}
	case utility::sql_database::DatabaseMode::sqlite3 : {
		statement_string = "INSERT OR IGNORE into "+table_name+"(";
		for ( size_t i=0; i<column_names.size(); i++ ) {
			statement_string+=column_names[i];
			if ( i != column_names.size()-1 ) {
				statement_string+=",";
			}
		}

		statement_string+=") VALUES(";
		for ( size_t i=0; i<values.size(); i++ ) {
			statement_string += "?";
			if ( i != column_names.size()-1 ) {
				statement_string+=",";
			}
		}
		statement_string+=");";

		statement stmt ( safely_prepare_statement(statement_string, db_session) );
		for ( size_t i=0; i<values.size(); i++ ) {
			stmt.bind(i+1, values[i]);
		}

		safely_write_to_database(stmt);
		break;
	}
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

void write_schema_to_database(
	string schema_str,
	sessionOP db_session)
{
	boost::char_separator< char > sep(";");
	boost::tokenizer< boost::char_separator< char > > tokens( schema_str, sep );
	BOOST_FOREACH ( std::string const & stmt_str, tokens ) {
		string trimmed_stmt_str(utility::trim(stmt_str, " \n\t"));
		if ( trimmed_stmt_str.size() ) {
			try{
				statement stmt = (*db_session) << trimmed_stmt_str + ";";
				safely_write_to_database(stmt);
			} catch (cppdb_error const & e) {
				TR.Error
					<< "ERROR reading schema \n"
					<< trimmed_stmt_str << std::endl;
				TR.Error << e.what() << std::endl;
				utility_exit();
			}
		}
	}
}


void
set_cache_size(
	sessionOP db_session,
	Size cache_size
) {

	if ( db_session->get_db_mode() == DatabaseMode::sqlite3 ) {
		stringstream stmt_ss;
		stmt_ss << "PRAGMA cache_size = " << cache_size << ";";
		statement stmt(safely_prepare_statement(stmt_ss.str(), db_session));
		safely_write_to_database(stmt);
	} else {
		TR
			<< "WARNING: Attempting to set database cache size "
			<< "for a database type for which this is currently not supported: "
			<< "'" << name_from_database_mode(db_session->get_db_mode()) << "'." << std::endl;
	}
}

std::string make_compound_statement(
	std::string const & table_name,
	std::vector<std::string> const & column_names,
	platform::Size const & row_count)
{
	std::string table_definition = table_name + " (" + utility::join(column_names,",") + ")";
	std::string value_list;
	platform::Size column_count = column_names.size();
	for ( platform::Size i = 0; i < row_count; ++i ) {
		std::string row_block= "(?";
		for ( platform::Size j = 1; j < column_count; ++j ) {
			row_block += ",?";
		}
		row_block += ")";

		value_list += row_block;
		if ( i != row_count-1 ) {
			value_list += ", ";
		}
	}

	return "INSERT INTO "+table_definition+ " VALUES " + value_list +";";
}

/// @detail build database connection from options in a tag, this is useful make sure the fields for
///constructing a database connection are consistent across different tags.
utility::sql_database::sessionOP
parse_database_connection(
	utility::tag::TagCOP tag
){
	using std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::inout;
	using namespace basic::options::OptionKeys;
	using utility::sql_database::DatabaseSessionManager;
	using namespace basic::resource_manager;

	bool separate_database = option[ inout::dbms::separate_db_per_mpi_process ]();
	int database_partition = option[ inout::dbms::database_partition]();

	if ( tag->hasOption("database_resource") ) {
		std::string database_resource = tag->getOption<string>("database_resource");
		if ( ! ResourceManager::get_instance()->has_resource_with_description( database_resource ) ) {
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a database_resource of '" + database_resource +
				"', but the ResourceManager doesn't have a resource with that description." );
		}
		return get_resource< utility::sql_database::session >( database_resource );
	}

	if ( tag->hasOption("database_resource_tag") ) {
		std::string database_resource_tag = tag->getOption<string>(
			"database_resource_tag");
		if ( ! ResourceManager::get_instance()->has_resource(
				database_resource_tag ) ) {
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a database_resource_tag of '" + database_resource_tag +
				"', but the ResourceManager doesn't have a resource with that tag." );
		}
		utility::sql_database::sessionOP db_session(utility::pointer::dynamic_pointer_cast< utility::sql_database::session > (
			ResourceManager::get_instance()->find_resource(database_resource_tag)));
		if ( !db_session ) {
			stringstream err_msg;
			err_msg
				<< "You specified a database_resource_tag of '" + database_resource_tag + "', while the ResourceManager does have a resource with that tag, it couldn't cast into a database session.";
			throw utility::excn::EXCN_Msg_Exception(err_msg.str());
		}
		return db_session;
	}

	utility::sql_database::TransactionMode::e transaction_mode = utility::sql_database::transaction_mode_from_name(
		tag->getOption<string>("transaction_mode", "standard"));

	Size chunk_size;
	switch(transaction_mode){
	case(utility::sql_database::TransactionMode::none) :
		if ( tag->hasOption("chunk_size") ) {
			TR << "WARNING: You must specify 'transaction_mode=chunk' ";
			TR << "to use the 'chunk size' tag." << endl;
		}
		chunk_size=0;
		break;
	case(utility::sql_database::TransactionMode::standard) :
		if ( tag->hasOption("chunk_size") ) {
			TR << "WARNING: You must specify 'transaction_mode=chunk' ";
			TR << "to use the 'chunk size' tag." << endl;
		}
		chunk_size=0;
		break;
	case(utility::sql_database::TransactionMode::chunk) :
		if ( !tag->hasOption("chunk_size") ) {
			utility_exit_with_message(
				"Must specify chunk_size if using the chunk transaction mode");
		}
		chunk_size=
			tag->getOption<Size>("chunk_size");
		break;
	default :
		utility_exit_with_message(
			"Unrecognized transaction mode: '" +
			name_from_transaction_mode(transaction_mode) + "'");
	}

	utility::sql_database::DatabaseMode::e database_mode;

	if ( tag->hasOption("database_mode") ) {
		database_mode = utility::sql_database::database_mode_from_name(
			tag->getOption<string>("database_mode"));
	} else {
		database_mode = utility::sql_database::database_mode_from_name(
			option[dbms::mode]);
	}

	std::string database_name;
	if ( tag->hasOption("database_name") ) {
		database_name = tag->getOption<string>("database_name");
	} else {
		database_name = option[dbms::database_name];
	}

	std::string database_dir;
	if ( option[ out::path::db ].user() ) {
		database_dir = option[ out::path::db ]().path();
	} else if ( option[ inout::dbms::path ].user() ) {
		database_dir = option[ inout::dbms::path ]().path();
	} else if ( option[ out::path::all ].user() ) {
		database_dir = option[ out::path::all ]().path();
	} else {
		database_dir = "";
	}
	std::string database_path = database_dir + database_name;

	// Parse pq_schema
	if ( tag->hasOption("database_pq_schema") && (database_mode != utility::sql_database::DatabaseMode::postgres) ) {
		TR << "WARNING: You must specify 'database_mode=postgres' ";
		TR << "to use the 'database_pq_schema' tag." << endl;
	}

	std::string database_pq_schema;
	if ( tag->hasOption("database_pq_schema") ) {
		database_pq_schema = tag->getOption<string>("database_pq_schema");
	} else {
		database_pq_schema = option[dbms::pq_schema];
	}


	// Check for invalid tags
	if ( database_mode != utility::sql_database::DatabaseMode::sqlite3 ) {
		if ( tag->hasOption("database_separate_db_per_mpi_process") ) {
			TR << "WARNING: You must specify 'database_mode=sqlite3' ";
			TR << "to use the 'database_separate_db_per_mpi_process' tag." << endl;
		}

		if ( tag->hasOption("database_partition") ) {
			TR << "WARNING: You must specify 'database_mode=sqlite3' ";
			TR << "to use the 'database_partition' tag." << endl;
		}

		if ( tag->hasOption("database_read_only") ) {
			TR << "WARNING: You must specify 'database_mode=sqlite3' ";
			TR << "to use the 'database_read_only tag." << endl;
		}
	}

	if ( database_mode == utility::sql_database::DatabaseMode::sqlite3 ) {
		if ( tag->hasOption("database_separate_db_per_mpi_process") && tag->hasOption("database_partition") ) {
			TR << "WARNING: 'database_separate_db_per_mpi_process' and 'database_partition' tags are mutually exclusive, using 'database_separate_db_per_mpi_process'." << endl;
		}

		if ( tag->hasOption("database_host") ) {
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_host' tag." << endl;
		}

		if ( tag->hasOption("database_user") ) {
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_user' tag." << endl;
		}

		if ( tag->hasOption("database_password") ) {
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_password' tag." << endl;
		}

		if ( tag->hasOption("database_port") ) {
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_port' tag." << endl;
		}
	}

	switch(database_mode){

	case utility::sql_database::DatabaseMode::sqlite3 :
		return DatabaseSessionManager::get_instance()->get_db_session(
			database_mode, transaction_mode, chunk_size,
			database_path, "", "", "", "", 0,
			tag->getOption("database_read_only", false),
			resolve_db_partition(
			tag->getOption("database_separate_db_per_mpi_process", separate_database),
			tag->getOption("database_partition", database_partition)));

	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres : {

		std::string database_host;
		if ( !tag->hasOption("database_host") ) {
			if ( !option[dbms::host].user() ) {
				utility_exit_with_message(
					"WARNING: To connect to a postgres or mysql database you must set"
					" the database_host tag or specify -dbms:host on the command line.");
			} else {
				database_host=option[dbms::host];
			}
		} else {
			database_host=tag->getOption<string>("database_host");
		}

		std::string database_user;
		if ( !tag->hasOption("database_user") ) {
			if ( !option[dbms::user].user() ) {
				utility_exit_with_message(
					"WARNING: To connect to a postgres or mysql database you must set"
					"the database_user tag or specify -dbms:user on the command line.");
			} else {
				database_user=option[dbms::user];
			}
		} else {
			database_user=tag->getOption<string>("database_user");
		}

		std::string database_password;
		if ( !tag->hasOption("database_password") ) {
			if ( !option[dbms::password].user() ) {
				utility_exit_with_message(
					"WARNING: To connect to a postgres or mysql database you must set"
					"the database_password tag or specify -dbms:password on the command line.");
			} else {
				database_password=option[dbms::password];
			}
		} else {
			database_password=tag->getOption<string>("database_password");
		}

		Size database_port;
		if ( !tag->hasOption("database_port") ) {
			if ( !option[dbms::port].user() ) {
				utility_exit_with_message(
					"WARNING: To connect to a postgres or mysql database you must set"
					"the database_port tag or specify -dbms:port on the command line.");
			} else {
				database_port=option[dbms::port];
			}
		} else {
			database_port=tag->getOption<Size>("database_port");
		}


		return DatabaseSessionManager::get_instance()->get_db_session(
			database_mode, transaction_mode, chunk_size,
			database_path,
			database_pq_schema,
			database_host,
			database_user,
			database_password,
			database_port);
	}
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(database_mode) + "'");
	}
	return 0;
}

}
}
