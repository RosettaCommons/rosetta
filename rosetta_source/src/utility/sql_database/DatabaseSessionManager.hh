// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/DatabaseSessionManager.hh
/// @author Matthew O'Meara
/// @author Sam Deluca
/// @author Chris Miles


#ifndef INCLUDE_utility_sql_database_DatabaseSessionManager_HH
#define INCLUDE_utility_sql_database_DatabaseSessionManager_HH

// Unit Headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <utility/pointer/ReferenceCount.hh>

// Boost Headers
#include <boost/scoped_ptr.hpp>

// C++ Headers
#include <string>

// External
#include <cppdb/frontend.h>

namespace utility {
namespace sql_database {

class session : public cppdb::session, public utility::pointer::ReferenceCount {};

class DatabaseSessionManager {

protected:

	// Private constructor to make it singleton managed
	DatabaseSessionManager();
	DatabaseSessionManager( const DatabaseSessionManager & src );

public:

	// Warning this is not called because of the singleton pattern
	virtual ~DatabaseSessionManager();

	///@brief return singleton instance of session manager
	static
	DatabaseSessionManager *
	get_instance();

	// for SQLite the database_url is the file path for the database
	sessionOP
	get_session(
		std::string const & db_fname);

	// overloaded get_session function for mysql
	sessionOP
	get_session(
		std::string const & host, std::string const & user, std::string const & password, std::string const & database,int const & port);

private:

#ifndef MULTITHREADED
	static boost::scoped_ptr< DatabaseSessionManager > instance_;
#endif
};



}
}



#endif
