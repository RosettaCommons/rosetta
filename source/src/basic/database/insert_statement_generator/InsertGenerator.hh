// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/insert_statement_generator/InsertGenerator.hh
///
/// @brief Insert statement generator
/// @author Sam DeLuca

#ifndef INCLUDED_basic_database_insert_statement_generator_InsertGenerator_HH
#define INCLUDED_basic_database_insert_statement_generator_InsertGenerator_HH

#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

#include <map>
#include <vector>


namespace basic {
namespace database {
namespace insert_statement_generator {

class InsertGenerator
{
public:
	InsertGenerator(std::string const & table_name);
	void add_column(std::string const & column_name);
	void add_row(std::vector<RowDataBaseOP> const & row);

	void write_to_database(utility::sql_database::sessionOP db_session);
	void write_to_database(utility::sql_database::sessionOP db_session, long long & last_sequence_id, std::string const & sequence_name);

private:
	void write_to_database_sequential(utility::sql_database::sessionOP db_session);
	void write_to_database_chunked(utility::sql_database::sessionOP db_session, platform::Size chunksize);

	void bind_row_data(cppdb::statement & statement, platform::Size row_start_index, platform::Size row_end_index);

private:
	std::string table_name_;

	std::vector<std::vector<RowDataBaseOP> > row_list_;
	std::vector<std::string> column_list_;

	std::map<std::string,platform::Size> column_index_map_;
};

}
}
}

#endif

