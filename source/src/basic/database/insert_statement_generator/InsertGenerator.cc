// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/insert_statement_generator/InsertGenerator.cc
///
/// @brief Insert statement generator
/// @author Sam DeLuca

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/sql_utils.hh>

#include <utility/string_util.hh>
#include <utility/sql_database/types.hh>
#include <cppdb/frontend.h>

namespace basic {
namespace database {
namespace insert_statement_generator {

InsertGenerator::InsertGenerator(std::string const & table_name) : table_name_(table_name)
{

}

void InsertGenerator::add_column(std::string const & column_name)
{
	column_list_.push_back(column_name);

	// Index into prepared statement for bind method, 1-based index.
	platform::Size column_index = column_list_.size();
	column_index_map_.insert(std::make_pair(column_name,column_index));
}

void InsertGenerator::add_row(std::vector<RowDataBaseOP> const & row)
{
	row_list_.push_back(row);
}

void
InsertGenerator::write_to_database(
	utility::sql_database::sessionOP db_session
) {
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		write_to_database_sequential(db_session);
		break;
	case utility::sql_database::DatabaseMode::mysql :
		write_to_database_chunked(db_session, 5000);
		break;
	case utility::sql_database::DatabaseMode::postgres :
		write_to_database_chunked(db_session, 300);
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

void
InsertGenerator::write_to_database(
	utility::sql_database::sessionOP db_session,
	long long & last_insert_id,
	std::string const & sequence_name
){
	if ( row_list_.size() > 1 ) {
		utility_exit_with_message("Unable to write to database with more than 1 row when requesting auto-increment sequence id.");
	}

	std::string statement_string = make_compound_statement(table_name_, column_list_, 1);
	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));

	bind_row_data(statement, 0, 1);
	basic::database::safely_write_to_database(statement);

	//Can't use cppdb::statement::last_insert_id() with postgres backend
	//last_insert_id = statement.last_insert_id();
	last_insert_id = statement.sequence_last(sequence_name);
}

void InsertGenerator::write_to_database_sequential(utility::sql_database::sessionOP db_session)
{
	std::string statement_string = make_compound_statement(table_name_, column_list_, 1);
	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( platform::Size i = 0; i < row_list_.size(); ++i ) {
		bind_row_data(statement, i, i+1);
		basic::database::safely_write_to_database(statement);
	}
}

void
InsertGenerator::write_to_database_chunked(
	utility::sql_database::sessionOP db_session,
	platform::Size chunk_size)
{
	std::string statement_string;
	cppdb::statement statement;

	platform::Size total_rows = row_list_.size();
	platform::Size row_start_index = 0;
	platform::Size row_end_index = 0;
	platform::Size remaining_rows = total_rows;
	//platform::Size column_count = column_list_.size();

	while ( remaining_rows > 0 )
			{
		if ( !statement.empty() ) {
			statement.reset();
		}

		platform::Size chunk = 0;
		if ( remaining_rows > chunk_size ) {
			chunk = chunk_size;
		} else {
			chunk = remaining_rows;
		}
		row_end_index = row_start_index+chunk;

		statement_string = basic::database::make_compound_statement(table_name_,column_list_,chunk);
		statement = basic::database::safely_prepare_statement(statement_string,db_session);

		bind_row_data(statement, row_start_index, row_end_index);

		row_start_index = row_end_index;
		basic::database::safely_write_to_database(statement);
		remaining_rows -= chunk;
	}
}

void InsertGenerator::bind_row_data(
	cppdb::statement & statement,
	platform::Size row_start_index,
	platform::Size row_end_index)
{
	platform::Size column_count = column_list_.size();

	for ( platform::Size i = row_start_index; i < row_end_index; ++i ) {
		std::vector<RowDataBaseOP> current_row = row_list_[i];
		for ( std::vector<RowDataBaseOP>::iterator column_it = current_row.begin(); column_it != current_row.end(); ++column_it ) {
			std::map<std::string,platform::Size>::const_iterator it(column_index_map_.find((*column_it)->get_column_name()));
			if ( it == column_index_map_.end() ) {
				utility_exit_with_message(table_name_ + " does not contain column " + (*column_it)->get_column_name() + " check for typos in your features reporter");
			}
			platform::Size base_column_index = it->second;
			platform::Size column_index = column_count*(i-row_start_index)+base_column_index;
			(*column_it)->bind_data(column_index,statement);
		}
	}
}

}
}
}
