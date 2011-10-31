// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/database/sql_utils.hh
/// @author Sam DeLuca

#ifndef INCLUDED_basic_database_sql_utils_HH_
#define INCLUDED_basic_database_sql_utils_HH_

#include <utility/sql_database/DatabaseSessionManager.hh>

namespace basic {
namespace database {

utility::sql_database::sessionOP get_db_session(
	std::string const & db_name,
	bool const readonly = false,
	bool const separate_db_per_mpi_process = false);

utility::sql_database::sessionOP get_db_session(
	std::string const & db_name,
	std::string const & db_mode,
	bool const readonly = false,
	bool const separate_db_per_mpi_process = false);

cppdb::statement safely_prepare_statement(std::string const & statement_string, utility::sql_database::sessionOP db_session);

void safely_write_to_database(cppdb::statement & statement);

cppdb::result safely_read_from_database(cppdb::statement & statement);


bool
table_exists(
	utility::sql_database::sessionOP db_session,
	std::string const & table_name);

}
}

#endif /* SQL_UTILS_HH_ */
