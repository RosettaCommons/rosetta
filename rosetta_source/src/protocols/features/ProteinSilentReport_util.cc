// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/ProteinSilentReport_util.cc
/// @author Sam DeLuca

#include <protocols/features/ProteinSilentReport_util.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>


namespace protocols {
namespace features {

core::Size get_current_structure_count_by_input_tag(
	utility::sql_database::sessionOP db_session,
	std::string const & input_tag)
{
	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	count(*)\n"
				"FROM\n"
				"	structures\n"
				"WHERE\n"
				"	structures.input_tag=?;" << input_tag;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
	core::Size count = 0;
	if(result.next())
	{
		result >> count;
	}
	return count;
}

core::Size get_score_type_id_from_score_term(
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	std::string const & score_term
)
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	score_type_id\n"
				"FROM\n"
				"	score_types\n"
				"WHERE\n"
				"	protocol_id=?\n"
				"AND\n"
				"	score_type_name=?;" << protocol_id << score_term;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
	core::Size score_type_id = 0;
	if(result.next())
	{
		result >> score_type_id;
	}
	return score_type_id;
}

core::Size get_struct_id_with_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag )
{
	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	job_string_real_data.struct_id\n"
				"FROM\n"
				"	job_string_real_data\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	job_string_real_data.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	job_string_real_data.data_key== ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	job_string_real_data.data_value\n"
				"LIMIT 1;" << score_term << input_tag;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Size get_struct_id_with_highest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag )
{
	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	job_string_real_data.struct_id\n"
				"FROM\n"
				"	job_string_real_data\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	job_string_real_data.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	job_string_real_data.data_key== ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	job_string_real_data.data_value DESC\n"
				"LIMIT 1;" << score_term << input_tag;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Size get_struct_id_with_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag )
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	structure_scores.struct_id\n"
				"FROM\n"
				"	structure_scores\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	structure_scores.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	structure_scores.score_type_id == ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	structure_scores.score_value\n"
				"LIMIT 1;" << score_type_id << input_tag;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Size get_struct_id_with_highest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag )
{
	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	structure_scores.struct_id\n"
				"FROM\n"
				"	structure_scores\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	structure_scores.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	structure_scores.score_type_id == ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	structure_scores.score_value DESC\n"
				"LIMIT 1;" << score_type_id << input_tag;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Size get_struct_id_with_nth_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag,
	core::Size const & cutoff_index)
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	job_string_real_data.struct_id\n"
				"FROM\n"
				"	job_string_real_data\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	job_string_real_data.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	job_string_real_data.data_key== ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	job_string_real_data.data_value\n"
				"LIMIT ?,1;" << score_term << input_tag << cutoff_index-1;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Size get_struct_id_with_nth_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag,
	core::Size const & cutoff_index)
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	structure_scores.struct_id\n"
				"FROM\n"
				"	structure_scores\n"
				"INNER JOIN\n"
				"	structures\n"
				"ON\n"
				"	structure_scores.struct_id == structures.struct_id\n"
				"WHERE\n"
				"	structure_scores.score_type_id == ?\n"
				"AND\n"
				"	structures.input_tag == ?\n"
				"ORDER BY\n"
				"	structure_scores.score_value\n"
				"LIMIT ?,1;" << score_type_id << input_tag << cutoff_index-1;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Size struct_id = 0;
	if(result.next())
	{
		result >> struct_id;
	}
	return struct_id;
}

core::Real get_score_for_struct_id_and_score_term_from_job_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	std::string const & score_term
	)
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	job_string_real_data.data_value\n"
				"FROM\n"
				"	job_string_real_data\n"
				"WHERE\n"
				"	job_string_real_data.data_key == ?\n"
				"AND\n"
				"	job_string_real_data.struct_id == ?\n;" << score_term <<struct_id;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Real score_value = 0.0;
	result >> score_value;
	return score_value;

}

core::Real get_score_for_struct_id_and_score_term_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	core::Size const & score_type_id)
{

	cppdb::result result;
	while(true)
	{
		try
		{
			result = (*db_session) <<
				"SELECT\n"
				"	structure_scores.score_value\n"
				"FROM\n"
				"	structure_scores\n"
				"WHERE\n"
				"	structure_scores.struct_id == ?\n"
				"AND\n"
				"	structure_scores.score_type_id == ?;" << struct_id << score_type_id;
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}

	core::Real score_value = 0.0;
	result >> score_value;
	return score_value;
}


}
}
