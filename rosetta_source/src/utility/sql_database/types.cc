// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/types.cc
/// @author Matthew O'Meara

#include <utility/sql_database/types.hh>
#include <utility/exit.hh>
#include <platform/types.hh>
#include <string>
#include <sstream>

namespace utility {
namespace sql_database {

using std::string;
using std::stringstream;

DatabaseMode::e
database_mode_from_name(
  std::string database_mode
) {
  if(!database_mode.compare("sqlite3")){
    return DatabaseMode::sqlite3;
  } else if(!database_mode.compare("mysql")){
    return DatabaseMode::mysql;
  } else if(!database_mode.compare("postgres")){
    return DatabaseMode::postgres;
  } else {
    utility_exit_with_message(
      "Unrecognized database mode: '" + database_mode + "'");
  }
  return DatabaseMode::sqlite3; // make compiler happy
}

std::string
name_from_database_mode(
  DatabaseMode::e database_mode
) {
  switch(database_mode){
  case DatabaseMode::sqlite3:
    return "sqlite3";
  case DatabaseMode::mysql:
    return "mysql";
  case DatabaseMode::postgres:
    return "postgres";
  default:
    stringstream err_msg;
    err_msg
      << "Unrecognized databse mode: '"
      << static_cast<platform::Size>(database_mode) << "'"; 
    utility_exit_with_message(err_msg.str());
  }
}

}
}
