// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
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

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace basic {
namespace database {
namespace insert_statement_generator{

InsertGenerator::InsertGenerator(std::string const & table_name) : table_name_(table_name)
{

}

void InsertGenerator::add_column(std::string const & column_name)
{
	platform::Size column_index = column_index_map_.size()+1;
	column_index_map_.insert(std::make_pair(column_name,column_index));
	index_column_map_.insert(std::make_pair(column_index,column_name));

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
	case utility::sql_database::DatabaseMode::sqlite3:
		write_to_database_sqlite(db_session);
		break;
	case utility::sql_database::DatabaseMode::mysql:
		write_to_database_mysql(db_session);
		break;
	case utility::sql_database::DatabaseMode::postgres:
		write_to_database_postgres(db_session);
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

void InsertGenerator::write_to_database_sqlite(utility::sql_database::sessionOP db_session)
{
	std::string columns = make_column_list();
	std::string placeholder_block =  "(?";
	for(platform::Size j = 2; j <= column_index_map_.size(); ++j)
	{
		placeholder_block += ",?";
	}
	placeholder_block += ")";

	std::string statement_string = "INSERT INTO "+table_name_+" "+columns + " VALUES "+placeholder_block+";";
	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
	for(std::vector<std::vector<RowDataBaseOP> >::iterator it = row_list_.begin(); it != row_list_.end();++it)
	{
		for(std::vector<RowDataBaseOP>::iterator column_it = it->begin(); column_it != it->end(); ++column_it)
		{
			std::map<std::string,platform::Size>::const_iterator it(column_index_map_.find((*column_it)->get_column_name()));
			if(it == column_index_map_.end())
			{
				utility_exit_with_message(table_name_ + " does not contain column " + (*column_it)->get_column_name() + " check for typos in your features reporter");
			}
			platform::Size index = it->second;
			(*column_it)->bind_data(index,statement);
		}
		basic::database::safely_write_to_database(statement);
	}

}


void
InsertGenerator::write_to_database_mysql(
	utility::sql_database::sessionOP db_session,
	platform::Size chunk_size
) {
	std::vector<std::string> column_names;
	for(platform::Size i = 1; i <= index_column_map_.size();++i)
	{
		column_names.push_back(index_column_map_[i]);
	}

	platform::Size total_rows = row_list_.size();
	platform::Size row_start_index = 0;
	platform::Size row_end_index = 0;
	platform::Size remaining_rows = total_rows;
	platform::Size column_count = column_index_map_.size();
	while(remaining_rows > 0)
	{
		platform::Size chunk = 0;
		if(remaining_rows > chunk_size)
		{
			chunk = chunk_size;
		}else
		{
			chunk = remaining_rows;
		}
		row_end_index = row_start_index+chunk;

		std::string statement_string(basic::database::make_compound_statement(table_name_,column_names,chunk));
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));

		for(platform::Size i = row_start_index; i < row_end_index;++i)
		{
			std::vector<RowDataBaseOP> current_row = row_list_[i];
			for(std::vector<RowDataBaseOP>::iterator column_it = current_row.begin(); column_it != current_row.end();++column_it)
			{
				std::map<std::string,platform::Size>::const_iterator it(column_index_map_.find((*column_it)->get_column_name()));
				if(it == column_index_map_.end())
				{
					utility_exit_with_message(table_name_ + " does not contain column " + (*column_it)->get_column_name() + " check for typos in your features reporter");
				}
				platform::Size base_column_index = it->second;
				platform::Size column_index = column_count*(i-row_start_index)+base_column_index;
				(*column_it)->bind_data(column_index,statement);
			}
		}
		row_start_index = row_end_index;
		basic::database::safely_write_to_database(statement);
		remaining_rows -= chunk;
	}
}

void
InsertGenerator::write_to_database_postgres(
	utility::sql_database::sessionOP db_session,
	platform::Size chunk_size
) {
	write_to_database_mysql(db_session, chunk_size);
}

std::string InsertGenerator::make_column_list() const
{
	std::string column_list = "(";
	for(platform::Size i = 1; i <= index_column_map_.size();++i)
	{
		std::string column = index_column_map_.find(i)->second;
		if(column_list.size() == 1)
		{
			column_list+= column;
		}else
		{
			column_list +=","+column;
		}
	}
	column_list+=")";

	return column_list;
}

}
}
}

