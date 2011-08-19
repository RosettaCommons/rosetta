// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/sqlite3_connection_manager.cc
/// @brief  Easy Access to sqlite3 database; manage connections
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifdef DB_SQLITE3

#ifndef INCLUDED_utility_sql_database_sqlite3_connection_manager_fwd_hh
#define INCLUDED_utility_sql_database_sqlite3_connection_manager_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace sql_database {

// Forward declarations
class Sqlite3ConnectionManager;
typedef utility::pointer::owning_ptr< Sqlite3ConnectionManager > Sqlite3ConnectionManagerOP;
typedef utility::pointer::owning_ptr< Sqlite3ConnectionManager const > Sqlite3ConnectionManagerCOP;


} // namespace sql_database
} // namespace utility

#endif // INCLUDED_utility_sql_database_sqlite3_connection_manager_FWD_HH

#endif // DB_SQLITE3
