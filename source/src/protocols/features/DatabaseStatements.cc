// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/features/ProteinSilentReport_util.cc
/// @author Sam DeLuca

#include <protocols/features/DatabaseStatements.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>

#include <sstream>

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
	if ( result.next() ) {
		result >> value;
		return value;
	} else {
		throw utility::excn::EXCN_Msg_Exception("No result found");
	}
}

cppdb::statement get_structure_count_statement(
	utility::sql_database::sessionOP db_session,
	std::string const & input_tag,
	core::Size const & protocol_id
){
	if ( input_tag.empty() ) {
		std::string statement_string =
			"SELECT\n"
			"\tcount(*)\n"
			"FROM\n"
			"\tstructures"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tbatches.protocol_id = ?;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,protocol_id);
		return statement;
	} else {
		std::string statement_string =
			"SELECT\n"
			"\tcount(*)\n"
			"FROM\n"
			"\tstructures\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tstructures.input_tag=?\n"
			"AND\n"
			"\tbatches.protocol_id=?;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,input_tag);
		statement.bind(2,protocol_id);
		return statement;
	}
}

core::Size get_score_type_id_from_score_term(
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	std::string const & score_term
)
{

	std::string statement_string =
		"SELECT\n"
		"score_type_id\n"
		"FROM\n"
		"score_types\n"
		"INNER JOIN\n"
		"batches ON score_types.batch_id = batches.batch_id\n"
		"WHERE\n"
		"batches.protocol_id=?\n"
		"AND\n"
		"score_type_name=?\n";

	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
	statement.bind(1,protocol_id);
	statement.bind(2,score_term);
	try{
		return get_something_from_database(statement, core::Size());
	}catch(utility::excn::EXCN_Msg_Exception &){
		throw;
		//utility_exit_with_message("No score_term "+score_term+" with protocol_id "+utility::to_string(protocol_id));
	}
}

cppdb::statement get_nth_lowest_score_from_job_data_statement(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	core::Size const & cutoff_index,
	std::string const & input_tag,
	core::Size const & protocol_id
){
	if ( input_tag.empty() ) {
		std::string statement_string =
			"SELECT\n"
			"\tjob_string_real_data.struct_id\n"
			"FROM\n"
			"\tjob_string_real_data\n"
			"INNER JOIN\n"
			"\tstructures ON job_string_real_data.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tjob_string_real_data.data_key = ? AND batches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tjob_string_real_data.data_value\n"
			"LIMIT ?,1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_term);
		statement.bind(2,protocol_id);
		statement.bind(3,cutoff_index-1);
		return statement;
	} else {

		std::string statement_string =
			"SELECT\n"
			"\tjob_string_real_data.struct_id\n"
			"FROM\n"
			"\tjob_string_real_data\n"
			"INNER JOIN\n"
			"\tstructures ON job_string_real_data.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tjob_string_real_data.data_key = ?\n"
			"AND\n"
			"\tstructures.input_tag = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tjob_string_real_data.data_value\n"
			"LIMIT ?,1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_term);
		statement.bind(2,input_tag);
		statement.bind(3,protocol_id);
		statement.bind(4,cutoff_index-1);
		return statement;
	}
}

cppdb::statement get_nth_lowest_score_from_score_data_statement(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	core::Size const & cutoff_index,
	std::string const & input_tag,
	core::Size const & protocol_id
){
	if ( input_tag.empty() ) {
		std::string statement_string =
			"SELECT\n"
			"\tstructure_scores.struct_id\n"
			"FROM\n"
			"\tstructure_scores\n"
			"INNER JOIN\n"
			"\tstructures ON structure_scores.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tstructure_scores.score_type_id = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tstructure_scores.score_value\n"
			"LIMIT ?,1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_type_id);
		statement.bind(2,protocol_id);
		statement.bind(3,cutoff_index-1);
		return statement;
	} else {
		std::string statement_string =
			"SELECT\n"
			"\tstructure_scores.struct_id\n"
			"FROM\n"
			"\tstructure_scores\n"
			"INNER JOIN\n"
			"\tstructures ON structure_scores.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tstructure_scores.score_type_id = ?\n"
			"AND\n"
			"\tstructures.input_tag = ?\n"
			"AND\n"
			"\tbatches.protocol_id =?\n"
			"ORDER BY\n"
			"\tstructure_scores.score_value\n"
			"LIMIT ?,1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_type_id);
		statement.bind(2,input_tag);
		statement.bind(3,protocol_id);
		statement.bind(4,cutoff_index-1);
		return statement;
	}
}

cppdb::statement get_highest_score_from_job_data_statement(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	std::string const & input_tag,
	core::Size const & protocol_id
){
	if ( input_tag.empty() ) {
		std::string statement_string =
			"SELECT\n"
			"\tjob_string_real_data.struct_id\n"
			"FROM\n"
			"\tjob_string_real_data\n"
			"INNER JOIN\n"
			"\tstructures ON job_string_real_data.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tjob_string_real_data.data_key = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tjob_string_real_data.data_value DESC\n"
			"LIMIT 1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_term);
		statement.bind(2,protocol_id);
		return statement;
	} else {
		std::string statement_string =
			"SELECT\n"
			"\tjob_string_real_data.struct_id\n"
			"FROM\n"
			"\tjob_string_real_data\n"
			"INNER JOIN\n"
			"\tstructures ON job_string_real_data.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tjob_string_real_data.data_key = ?\n"
			"AND\n"
			"\tstructures.input_tag = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tjob_string_real_data.data_value DESC\n"
			"LIMIT 1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_term);
		statement.bind(2,input_tag);
		statement.bind(3,protocol_id);
		return statement;
	}
}

cppdb::statement get_highest_score_from_score_data_statement(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	std::string const & input_tag,
	core::Size const & protocol_id
){
	if ( input_tag.empty() ) {
		std::string statement_string =
			"SELECT\n"
			"\tstructure_scores.struct_id\n"
			"FROM\n"
			"\tstructure_scores\n"
			"INNER JOIN\n"
			"\tstructures ON structure_scores.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tstructure_scores.score_type_id = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tstructure_scores.score_value DESC\n"
			"LIMIT 1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_type_id);
		statement.bind(2,protocol_id);
		return statement;
	} else {
		std::string statement_string =
			"SELECT\n"
			"\tstructure_scores.struct_id\n"
			"FROM\n"
			"\tstructure_scores\n"
			"INNER JOIN\n"
			"\tstructures ON structure_scores.struct_id = structures.struct_id\n"
			"INNER JOIN\n"
			"\tbatches ON structures.batch_id = batches.batch_id\n"
			"WHERE\n"
			"\tstructure_scores.score_type_id = ?\n"
			"AND\n"
			"\tstructures.input_tag = ?\n"
			"AND\n"
			"\tbatches.protocol_id = ?\n"
			"ORDER BY\n"
			"\tstructure_scores.score_value DESC\n"
			"LIMIT 1;";
		cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
		statement.bind(1,score_type_id);
		statement.bind(2,input_tag);
		statement.bind(3,protocol_id);
		return statement;
	}
}

////// Functions declared in the .hh file //////

core::Size get_current_structure_count(
	utility::sql_database::sessionOP db_session,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	cppdb::statement statement = get_structure_count_statement(db_session, input_tag,protocol_id);
	return get_something_from_database(statement, core::Size());
}

StructureID
get_struct_id_with_nth_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	core::Size const & cutoff_index,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	cppdb::statement statement = get_nth_lowest_score_from_job_data_statement(db_session, score_term, cutoff_index, input_tag,protocol_id);

	try{
		return get_something_from_database(statement, StructureID());
	}catch(utility::excn::EXCN_Msg_Exception){
		std::stringstream message;
		message << "No nth lowest " << score_term << " with input_tag " << input_tag << ", where n=" <<utility::to_string(cutoff_index);

		utility_exit_with_message(message.str());
	}

// shouldn't get here
	return StructureID();
}

StructureID
get_struct_id_with_nth_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	core::Size const & cutoff_index,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	cppdb::statement statement = get_nth_lowest_score_from_score_data_statement(db_session, score_type_id, cutoff_index, input_tag,protocol_id);

	try{
		return get_something_from_database(statement, StructureID());
	}catch(utility::excn::EXCN_Msg_Exception){
		utility_exit_with_message("No nth lowest score_type_id: "+utility::to_string(score_type_id)+", with input_tag "+input_tag+ ", where n="+utility::to_string(cutoff_index));
	}

// shouldn't get here
	return StructureID();
}

StructureID get_struct_id_with_lowest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	return get_struct_id_with_nth_lowest_score_from_job_data(db_session, score_term, 1, protocol_id, input_tag);
}

StructureID get_struct_id_with_lowest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	return get_struct_id_with_nth_lowest_score_from_score_data(db_session, score_type_id, 1, protocol_id, input_tag);
}

StructureID get_struct_id_with_highest_score_from_job_data(
	utility::sql_database::sessionOP db_session,
	std::string const & score_term,
	core::Size const & protocol_id,
	std::string const & input_tag )
{
	cppdb::statement statement = get_highest_score_from_job_data_statement(db_session, score_term, input_tag, protocol_id);

	try{
		return get_something_from_database(statement, StructureID());
	}catch(utility::excn::EXCN_Msg_Exception){
		utility_exit_with_message("No "+score_term+" with input_tag "+input_tag);
	}

// shouldn't get here
	return StructureID();
}

StructureID get_struct_id_with_highest_score_from_score_data(
	utility::sql_database::sessionOP db_session,
	core::Size const & score_type_id,
	core::Size const & protocol_id,
	std::string const & input_tag
){
	cppdb::statement statement = get_highest_score_from_score_data_statement(db_session, score_type_id, input_tag,protocol_id);

	try{
		return get_something_from_database(statement, StructureID());
	}catch(utility::excn::EXCN_Msg_Exception){
		utility_exit_with_message("No score_type_id: "+utility::to_string(score_type_id)+", with input_tag "+input_tag);
	}

// shouldn't get here
	return StructureID();
}

core::Real
get_score_for_struct_id_and_score_term_from_job_data(
	utility::sql_database::sessionOP db_session,
	StructureID const & struct_id,
	std::string const & score_term
){

	std::string statement_string =
		"SELECT\n"
		"\tjob_string_real_data.data_value\n"
		"FROM\n"
		"\tjob_string_real_data\n"
		"WHERE\n"
		"\tjob_string_real_data.data_key = ?\n"
		"AND\n"
		"\tjob_string_real_data.struct_id = ?\n;";
	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
	statement.bind(1,score_term);
	statement.bind(2,struct_id);

	try{
		return get_something_from_database(statement, core::Real());
	}catch(utility::excn::EXCN_Msg_Exception){
		std::stringstream message;
		message << "No " << score_term << " with struct_id " << struct_id;

		utility_exit_with_message(message.str());
	}

// shouldn't get here
	return 0;
}

core::Real
get_score_for_struct_id_and_score_term_from_score_data(
	utility::sql_database::sessionOP db_session,
	StructureID const & struct_id,
	core::Size const & score_type_id)
{

	std::string statement_string =
		"SELECT\n"
		"\tstructure_scores.score_value\n"
		"FROM\n"
		"\tstructure_scores\n"
		"WHERE\n"
		"\tstructure_scores.struct_id = ?\n"
		"AND\n"
		"\tstructure_scores.score_type_id = ?;";
	cppdb::statement statement(basic::database::safely_prepare_statement(statement_string,db_session));
	statement.bind(1,struct_id);
	statement.bind(2,score_type_id);

	try{
		return get_something_from_database(statement, core::Real());
	}catch(utility::excn::EXCN_Msg_Exception){
		std::stringstream message;
		message << "No score_type_id:" << score_type_id << ", with struct_id " << struct_id;
		utility_exit_with_message(message.str());
	}

// shouldn't get here
	return 0;
}

}
}
