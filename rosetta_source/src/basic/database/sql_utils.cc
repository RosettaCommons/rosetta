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

using std::string;
using utility::sql_database::sessionOP;
using utility::sql_database::DatabaseSessionManager;

namespace basic {
namespace database {

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

}
}
