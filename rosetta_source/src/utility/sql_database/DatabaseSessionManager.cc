// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/DatabaseSessionManager.cc
/// @author Matthew O'Meara
/// @author Sam Deluca
/// @author Chris Miles

#ifdef USEMPI
#include <mpi.h>
#endif

// Unit Headers
#include <utility/sql_database/DatabaseSessionManager.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>

// Boost Headers
#ifdef MULTITHREADED
#include <boost/thread/tss.hpp>
#else
#include <boost/scoped_ptr.hpp>
#endif

// C++ Headers
#include <string>
#include <sstream>

namespace utility {
namespace sql_database {

using utility::file::FileName;
using std::string;
using std::stringstream;
using cppdb::cppdb_error;

#ifdef MULTITHREADED
boost::thread_specific_pointer< DatabaseSessionManager > DatabaseSessionManager::instance_;
#else
boost::scoped_ptr< DatabaseSessionManager > DatabaseSessionManager::instance_;
#endif

DatabaseSessionManager *
DatabaseSessionManager::get_instance(){
    if( instance_.get() == 0 ){
        instance_.reset( new DatabaseSessionManager() );
    }
    return instance_.get();
}

DatabaseSessionManager::DatabaseSessionManager() {}

DatabaseSessionManager::DatabaseSessionManager(
const DatabaseSessionManager &
) {}

DatabaseSessionManager::~DatabaseSessionManager() {}

///details@ Currently only supports SQLite3 and mysqldatabase databases
/// if the separate_db_per_mpi_process appends "_<mpi_rank>" to the end of the database filename
/// This is useful when writing to an sqlite database not through the job distributor where locking causes problems
sessionOP
DatabaseSessionManager::get_session(
  std::string const & db_fname,
  bool const readonly /* = false */,
    bool const separate_db_per_mpi_process /* = false */
){
    sessionOP s(new session());

    try {
        string use_db_fname;

#ifdef USEMPI
        if(separate_db_per_mpi_process){
            int mpi_rank(0);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            stringstream buf; buf << FileName(db_fname).name() << "_" << mpi_rank;
            use_db_fname = buf.str();
        } else {
            use_db_fname = FileName(db_fname).name();
        }
#else
        use_db_fname = FileName(db_fname).name();
#endif


        if(readonly){
            s->open("sqlite3:mode=readonly;db="+use_db_fname);
        } else {
            s->open("sqlite3:db="+use_db_fname);
        }
    } catch (cppdb_error & e){
        std::stringstream error_msg;
        error_msg
            << "Failed to open database file '" << db_fname << "'"
            << (readonly ? " in readonly mode:" : ":") << std::endl
            << "\t" << e.what();
        utility_exit_with_message(error_msg.str());
    }
    return s;
}

sessionOP
DatabaseSessionManager::get_session(
    std::string const & db_mode,                                
    std::string const & host,
    std::string const & user,
    std::string const & password,
    std::string const & database,
    int const & port    
){

    sessionOP s(new session());
    std::string port_string(utility::to_string<int>(port));
    
  try {
       
      if(db_mode == "postgres"){
#ifndef USEPOSTGRES
          utility_exit_with_message("If you want to use a postgres database, build with extras=postgres");
#endif
          s->open("postgresql:user="+user+";dbname="+database+";port="+port_string);
      }
      else if(db_mode == "mysql"){          
#ifndef USEMYSQL
          utility_exit_with_message("If you want to use a mysql database, build with extras=mysql");
#endif
          s->open("mysql:host="+host+";user="+user+";password="+password+";database="+database+";port="+port_string+";opt_reconnect=1");
      }
  } catch (cppdb_error & e){
      std::stringstream error_msg;
      error_msg
      << "Failed to open database file '" << database << "'"
      <<  std::endl
      << "\t" << e.what();
      utility_exit_with_message(error_msg.str());
  }      
   
    return s;
//#else
//	  	utility_exit_with_message("You shouldn't be here, also if you want to use mysql specify extras=mysql when you build");
//		return NULL; // need to return something to compile on Windows VC++
//#endif
}

} // namespace
} // namespace

