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
using std::stringstream;

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

  // Currently only supports SQLite3 and mysqldatabase databases
  sessionOP
  DatabaseSessionManager::get_session(
	  std::string const & db_fname,
	  bool const readonly
	){
		sessionOP s(new session());
		if(readonly){
			s->open("sqlite3:mode=readonly;db="+FileName(db_fname).name());
		} else {
			// for non-readonly databases, each node writes connects to it's own database file indexed by rank
#ifdef USEMPI
			stringstream part_db_fname;
			int mpi_rank(0);
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

			part_db_fname << FileName(db_fname).name() << "_" << mpi_rank;
			s->open("sqlite3:db="+part_db_fname.str());
#else
			s->open("sqlite3:db="+FileName(db_fname).name());
#endif
		}
		return s;
  }

  sessionOP
  DatabaseSessionManager::get_session(
		std::string const & host,
		std::string const & user,
		std::string const & password,
		std::string const & database,
		int const & port
	){
#ifdef USEMYSQL
	  	sessionOP s(new session());
	  	std::string port_string(utility::to_string<int>(port));
	  	s->open("mysql:host="+host+";user="+user+";password="+password+";database="+database+";port="+port_string);
	  	return s;
#else
	  	utility_exit_with_message("You shouldn't be here, also if you want to use mysql specify extras=mysql when you build");
		return NULL; // need to return something to compile on Windows VC++
#endif
  }

} // namespace
} // namespace

