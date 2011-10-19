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

#include <protocols/features/DatabaseStatements.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>


namespace protocols {
namespace features {

////// Functions only available to this CC file //////

template <class T>
T get_something_from_database(
		cppdb::statement statement,
		T // dummy variable. I hate dummy variables, but that's the best I can do for now
){
	cppdb::result result(basic::database::safely_read_from_database(statement));
	T value;
	if(result.next())
	{
		result >> value;
	}
	return value;
}

cppdb::statement get_structure_count_statement(
		utility::sql_database::sessionOP db_session,
		std::string const & input_tag
){
	if( input_tag.empty()){
		return (*db_session) <<
				"SELECT\n"
				"	count(*)\n"
				"FROM\n"
				"	structures;";
	}else{
		return 	(*db_session) << "SELECT\n"
				"	count(*)\n"
				"FROM\n"
				"	structures\n"
				"WHERE\n"
				"	structures.input_tag=?;" << input_tag;
	}
}

core::Size get_score_type_id_from_score_term(
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	std::string const & score_term
)
{

	cppdb::statement statement = (*db_session) <<
		"SELECT\n"
		"	score_type_id\n"
		"FROM\n"
		"	score_types\n"
		"WHERE\n"
		"	protocol_id=?\n"
		"AND\n"
		"	score_type_name=?;" << protocol_id << score_term;

	return get_something_from_database(statement, core::Size());
}

cppdb::statement get_nth_lowest_score_from_job_data_statement(
		utility::sql_database::sessionOP db_session,
		std::string const & score_term,
		core::Size const & cutoff_index,
		std::string const & input_tag
){
	if(input_tag.empty()){
		return 	(*db_session) <<
				"SELECT\n"
				"	job_string_real_data.struct_id\n"
				"FROM\n"
				"	job_string_real_data\n"
				"WHERE\n"
				"	job_string_real_data.data_key== ?\n"
				"ORDER BY\n"
				"	job_string_real_data.data_value\n"
				"LIMIT ?,1;" << score_term << cutoff_index-1;
	}
	else{
		return (*db_session) <<
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
	}
}

cppdb::statement get_nth_lowest_score_from_score_data_statement(
		utility::sql_database::sessionOP db_session,
		core::Size const & score_type_id,
		core::Size const & cutoff_index,
		std::string const & input_tag
){
	if(input_tag.empty()){
		return (*db_session) <<
				"SELECT\n"
				"	structure_scores.struct_id\n"
				"FROM\n"
				"	structure_scores\n"
				"WHERE\n"
				"	structure_scores.score_type_id == ?\n"
				"ORDER BY\n"
				"	structure_scores.score_value\n"
				"LIMIT ?,1;" << score_type_id << cutoff_index-1;
	}
	else{
		return (*db_session) <<
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
	}
}

cppdb::statement get_highest_score_from_job_data_statement(
		utility::sql_database::sessionOP db_session,
		std::string const & score_term,
		std::string const & input_tag
){
	if(input_tag.empty()){
		return 	(*db_session) <<
				"SELECT\n"
				"	job_string_real_data.struct_id\n"
				"FROM\n"
				"	job_string_real_data\n"
				"WHERE\n"
				"	job_string_real_data.data_key== ?\n"
				"ORDER BY\n"
				"	job_string_real_data.data_value DESC\n"
				"LIMIT 1;" << score_term;
	}
	else{
		return (*db_session) <<
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
	}
}

cppdb::statement get_highest_score_from_score_data_statement(
		utility::sql_database::sessionOP db_session,
		core::Size const & score_type_id,
		std::string const & input_tag
){
	if(input_tag.empty()){
		return (*db_session) <<
				"SELECT\n"
				"	structure_scores.struct_id\n"
				"FROM\n"
				"	structure_scores\n"
				"ORDER BY\n"
				"	structure_scores.score_value DESC\n"
				"LIMIT 1;" << score_type_id;
	}
	else{
		return (*db_session) <<
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
	}
}

////// Functions declared in the .hh file //////

core::Size get_current_structure_count(
	utility::sql_database::sessionOP db_session,
	std::string const & input_tag
){
	cppdb::statement statement = get_structure_count_statement(db_session, input_tag);
	return get_something_from_database(statement, core::Size());
}

core::Size get_struct_id_with_nth_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	core::Size const & cutoff_index,
	std::string const & input_tag
){
	cppdb::statement statement = get_nth_lowest_score_from_job_data_statement(db_session, score_term, cutoff_index, input_tag);
	return get_something_from_database(statement, core::Size());
}

core::Size get_struct_id_with_nth_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	core::Size const & cutoff_index,
	std::string const & input_tag
){
	cppdb::statement statement = get_nth_lowest_score_from_score_data_statement(db_session, score_type_id, cutoff_index, input_tag);
	return get_something_from_database(statement, core::Size());
}

core::Size get_struct_id_with_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag
){
	return get_struct_id_with_nth_lowest_score_from_job_data(db_session, score_term, 1, input_tag);
}

core::Size get_struct_id_with_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag
){
	return get_struct_id_with_nth_lowest_score_from_score_data(db_session, score_type_id, 1, input_tag);
}

core::Size get_struct_id_with_highest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag )
{
	cppdb::statement statement = get_highest_score_from_job_data_statement(db_session, score_term, input_tag);
	return get_something_from_database(statement, core::Size());
}

core::Size get_struct_id_with_highest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag
){
	cppdb::statement statement = get_highest_score_from_score_data_statement(db_session, score_type_id, input_tag);
	return get_something_from_database(statement, core::Size());
}

core::Real get_score_for_struct_id_and_score_term_from_job_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	std::string const & score_term
){

	cppdb::statement statement = (*db_session) <<
		"SELECT\n"
		"	job_string_real_data.data_value\n"
		"FROM\n"
		"	job_string_real_data\n"
		"WHERE\n"
		"	job_string_real_data.data_key == ?\n"
		"AND\n"
		"	job_string_real_data.struct_id == ?\n;" << score_term <<struct_id;

	return get_something_from_database(statement, core::Real());
}

core::Real get_score_for_struct_id_and_score_term_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & struct_id,
	core::Size const & score_type_id)
{

	cppdb::statement statement = (*db_session) <<
		"SELECT\n"
		"	structure_scores.score_value\n"
		"FROM\n"
		"	structure_scores\n"
		"WHERE\n"
		"	structure_scores.struct_id == ?\n"
		"AND\n"
		"	structure_scores.score_type_id == ?;" << struct_id << score_type_id;

	return get_something_from_database(statement, core::Real());
}

}
}
