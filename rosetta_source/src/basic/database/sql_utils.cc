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
namespace basic {
namespace database {

utility::sql_database::sessionOP get_db_session(std::string const & db_name)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string db_mode(option[inout::database_mode]);
	return get_db_session(db_name,db_mode);
}

utility::sql_database::sessionOP get_db_session(std::string const & db_name, std::string const & db_mode)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//std::string db_mode(option[inout::database_mode]);
	if(db_mode == "sqlite3"){

		if(option[mysql::host].user() || option[mysql::user].user() || option[mysql::password].user() || option[mysql::port].user())
		{
			utility_exit_with_message("you have specified mysql server options, but are running in sqlite3 mode, spefcify -inout:database_mode mysql to rectify this.");
		}
		utility::sql_database::sessionOP db_session(utility::sql_database::DatabaseSessionManager::get_instance()->get_session(db_name));
		return db_session;
	}else if(db_mode == "mysql")
	{
#ifndef USEMYSQL
		utility_exit_with_message("if you want to use a mysql database, build with extras=mysql");
#endif
		if(option[mysql::host].user() && option[mysql::user].user() && option[mysql::password].user() && option[mysql::port].user())
		{
			std::string host(option[mysql::host]);
			std::string user(option[mysql::user]);
			std::string password(option[mysql::password]);
			platform::Size port(option[mysql::port]);


			utility::sql_database::sessionOP db_session(
				utility::sql_database::DatabaseSessionManager::get_instance()->get_session(host,user,password,db_name,port));
			return db_session;
		}else
		{
			utility_exit_with_message("yiou must specify the following options to use a mysql database: -mysql:host -mysql:user -mysql:password -mysql:port");
		}
	}else
	{
		utility_exit_with_message("you need to specify either 'mysql' or 'sqlite3' as a mode with -inout:database_mode.  You specified: "+db_mode);
	}
}

}
}
