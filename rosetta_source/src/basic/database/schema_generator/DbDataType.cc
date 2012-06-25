// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/DbDataType.cc
///
/// @brief DbDataType class for the schema generator framework
/// @author Tim Jacobs

#include <basic/database/schema_generator/DbDataType.hh>

//Basic
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Utility
#include <utility/exit.hh>
#include <utility/string_util.hh>

//C++
#include <string>

namespace basic{
namespace database{
namespace schema_generator{

/*********Base class*********/
DbDataType::DbDataType() {
	database_mode_ =
		utility::sql_database::database_mode_from_name(
			basic::options::option[basic::options::OptionKeys::inout::dbms::mode]);
}

std::string DbDataType::print() const{
	return type_string_;
}

/*********Text based data types*********/
DbText::DbText():
DbDataType()
{
	type_string_ = "TEXT";
}

DbText::DbText(int size):
DbDataType()
{

	std::string size_string = utility::to_string(size);
	switch(database_mode_){
	case utility::sql_database::DatabaseMode::sqlite3:
		type_string_ = "TEXT";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:
		type_string_ = "VARCHAR(" + size_string + ")";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(database_mode_) + "'");
	}

}

//Needed because MYSQL doesn't support variable length primary keys
DbTextKey::DbTextKey() :
	DbDataType()
{
	switch(database_mode_){
	case utility::sql_database::DatabaseMode::sqlite3:
		type_string_ = "TEXT";
		break;
	case utility::sql_database::DatabaseMode::postgres:
		type_string_ = "TEXT";
		break;
	case utility::sql_database::DatabaseMode::mysql:
		type_string_ = "VARCHAR(255)";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(database_mode_) + "'");
	}

}

/*********Integer data types*********/
DbInteger::DbInteger() :
	DbDataType()
{
	type_string_ = "INTEGER";
}

DbInteger::DbInteger(int size) :
	DbDataType()
{
	std::string size_string = utility::to_string(size);
	type_string_ = "INTEGER";
}

DbBigInt::DbBigInt() :
	DbDataType()
{

	switch(database_mode_){
	case utility::sql_database::DatabaseMode::sqlite3:
		type_string_ = "INTEGER";
		break;
	case utility::sql_database::DatabaseMode::postgres:
		type_string_ = "BIGINT";
		break;
	case utility::sql_database::DatabaseMode::mysql:
		type_string_ = "BIGINT";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(database_mode_) + "'");
	}
}

/// DOES NOT WORK WITH CPPDB - USE INTEGER INSTEAD
///*********Boolean data types*********/
//DbBoolean::DbBoolean():
//DbDataType()
//{
//	type_string_ = "BOOLEAN";
//}

DbReal::DbReal():
DbDataType()
{
	type_string_ = "REAL";
}

DbDouble::DbDouble():
DbDataType()
{

	switch(database_mode_){
	case utility::sql_database::DatabaseMode::sqlite3:
		type_string_ = "REAL";
		break;
	case utility::sql_database::DatabaseMode::postgres:
		type_string_ = "DOUBLE PRECISION";
		break;
	case utility::sql_database::DatabaseMode::mysql:
		//I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
		type_string_ = "DOUBLE";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(database_mode_) + "'");
	}
}
/*********UUID data types*********/
DbUUID::DbUUID() :
	DbDataType()
{

	switch(database_mode_){
	case utility::sql_database::DatabaseMode::sqlite3:
		type_string_ = "BLOB";
		break;
	case utility::sql_database::DatabaseMode::postgres:
		type_string_ = "UUID";
		break;
	case utility::sql_database::DatabaseMode::mysql:
		//I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
		type_string_ = "BINARY(16)";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(database_mode_) + "'");
	}
}

} // schema_generator
} // namespace database
} // namespace utility
