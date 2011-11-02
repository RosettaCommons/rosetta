// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/database/sql_utils.cc
/// @author Sam DeLuca

#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mysql.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/Tracer.hh>

#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/scoped_ptr.hpp>
#include <cppdb/frontend.h>


using std::string;
using utility::sql_database::sessionOP;
using utility::sql_database::DatabaseSessionManager;

namespace basic {
namespace database {

static basic::Tracer TR( "basic.database.sql_utils" );

utility::sql_database::sessionOP get_db_session(
	string const & db_name,
	bool const readonly /* = false */,
  bool const separate_db_per_mpi_process /* = false */
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	string db_mode(option[inout::database_mode]);
	return get_db_session(db_name, db_mode, readonly, separate_db_per_mpi_process);
}

sessionOP get_db_session(
	string const & db_name,
	string const & db_mode,
	bool const readonly /* = false */,
	bool const separate_db_per_mpi_process /* = false */
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//string db_mode(option[inout::database_mode]);
	if(db_mode == "sqlite3"){

		if(option[mysql::host].user() || option[mysql::user].user() || option[mysql::password].user() || option[mysql::port].user())
		{
			utility_exit_with_message("you have specified mysql server options, but are running in sqlite3 mode, specify -inout:database_mode mysql to rectify this.");
		}
		sessionOP db_session(DatabaseSessionManager::get_instance()->get_session(db_name, readonly, separate_db_per_mpi_process));
		return db_session;
	}else if(db_mode == "mysql")
	{
#ifndef USEMYSQL
		utility_exit_with_message("If you want to use a mysql database, build with extras=mysql");
#endif
		if(readonly){
			utility_exit_with_message("Restricting access to a mysql database is done at the user level rather that the connection level. So requesting a readonly connection cannot fullfilled.");
		}

		if(option[mysql::host].user() && option[mysql::user].user() && option[mysql::password].user() && option[mysql::port].user())
		{
			string host(option[mysql::host]);
			string user(option[mysql::user]);
			string password(option[mysql::password]);
			platform::Size port(option[mysql::port]);


			sessionOP db_session(
				DatabaseSessionManager::get_instance()->get_session(host,user,password,db_name,port));
			return db_session;
		}else
		{
			utility_exit_with_message("You must specify the following options to use a mysql database: -mysql:host -mysql:user -mysql:password -mysql:port");
		}
	}else
	{
		utility_exit_with_message("You need to specify either 'mysql' or 'sqlite3' as a mode with -inout:database_mode.  You specified: "+db_mode);
	}
}


cppdb::statement safely_prepare_statement(std::string const & statement_string, utility::sql_database::sessionOP db_session)
{
	cppdb::statement stmt;
	try
	{
		stmt = db_session->prepare(statement_string);
		return stmt;
	}catch(cppdb::cppdb_error error)
	{
		utility_exit_with_message(error.what());
	}
	return stmt; //there's no way this should happen
}

void safely_write_to_database(cppdb::statement & statement)
{
	while(true)
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
		}catch(cppdb::null_value_fetch & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::cppdb_error & except)
		{
			utility_exit_with_message(except.what());
			//TR <<except.what() <<std::endl;
			#ifndef WIN32
				usleep(10);
			#endif
			continue;
		}
	}
}

cppdb::result safely_read_from_database(cppdb::statement & statement)
{
	while(true)
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
			utility_exit_with_message(except.what());
			//TR <<except.what() <<std::endl;
			#ifndef WIN32
				usleep(10);
			#endif
			continue;
		}
	}
}




bool
table_exists(
	sessionOP db_session,
	string const & table_name
) {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	std::string db_name(basic::options::option[basic::options::OptionKeys::inout::database_filename].value_string());

	cppdb::statement stmt;
	if(db_mode == "sqlite3")
	{
		std::string statement_string = "SELECT name FROM sqlite_master WHERE name=?;";
		stmt = safely_prepare_statement(statement_string,db_session);
	}else if(db_mode == "mysql")
	{
		std::string statement_string = "SHOW TABLES WHERE Tables_in_"+db_name+" = ?;";
		stmt = safely_prepare_statement(statement_string,db_session);
	}else
	{
		utility_exit_with_message("unknown database mode");
	}

	stmt.bind(1,table_name);
	cppdb::result res = stmt.query();

	if(res.next()){
		return true;
	} else {
		return false;
	}
}

}
}
