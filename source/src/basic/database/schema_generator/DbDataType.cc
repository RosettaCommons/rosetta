// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/schema_generator/DbDataType.cc
///
/// @brief DbDataType class for the schema generator framework
/// @author Tim Jacobs

#include <basic/database/schema_generator/DbDataType.hh>
#include <memory>                                              // for shared...

//Utility
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//C++
#include <string>

namespace basic {
namespace database {
namespace schema_generator {

/*********Base class*********/
DbDataType::DbDataType():
	size_(0)
{}

DbDataType::DbDataType(int size):
	size_(size)
{}

/*********Text based data types*********/
DbText::DbText():
	DbDataType()
{}

DbText::DbText(int size):
	DbDataType(size)
{}

std::string
DbText::print(
	utility::sql_database::sessionOP db_session
) const {

	std::string size_string = utility::to_string(size_);
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "TEXT";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		if ( size_>0 ) {
			return "VARCHAR(" + size_string + ")";
		} else {
			return "TEXT";
		}
		break;
	case utility::sql_database::DatabaseMode::postgres :
		if ( size_>0 ) {
			return "VARCHAR(" + size_string + ")";
		} else {
			return "VARCHAR";
		}
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
	//appease the compiler
	return "";
}

//Needed because MYSQL doesn't support variable length primary keys
DbTextKey::DbTextKey() :
	DbDataType()
{}

std::string
DbTextKey::print(
	utility::sql_database::sessionOP db_session
) const {
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "TEXT";
		break;
	case utility::sql_database::DatabaseMode::postgres :
		return "TEXT";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		return "VARCHAR(255)";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

/*********Integer data types*********/
DbInteger::DbInteger() :
	DbDataType()
{}

std::string
DbInteger::print(
	utility::sql_database::sessionOP /*db_session*/
) const {
	return "INTEGER";
}

DbBigInt::DbBigInt() :
	DbDataType()
{}

std::string
DbBigInt::print(
	utility::sql_database::sessionOP db_session
) const {
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "INTEGER";
		break;
	case utility::sql_database::DatabaseMode::postgres :
		return "BIGINT";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		return "BIGINT";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

/// DOES NOT WORK WITH CPPDB - USE INTEGER INSTEAD
///*********Boolean data types*********/
//DbBoolean::DbBoolean():
//DbDataType()
//{
// return "BOOLEAN";
//}

DbReal::DbReal():
	DbDataType()
{}

std::string
DbReal::print(
	utility::sql_database::sessionOP db_session
) const {
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "REAL";
		break;
	case utility::sql_database::DatabaseMode::postgres :
		return "DOUBLE PRECISION";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		//I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
		return "REAL";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

DbDouble::DbDouble():
	DbDataType()
{}

std::string
DbDouble::print(
	utility::sql_database::sessionOP db_session
) const {

	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "REAL";
		break;
	case utility::sql_database::DatabaseMode::postgres :
		return "DOUBLE PRECISION";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		//I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
		return "DOUBLE";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}
/*********UUID data types*********/
DbUUID::DbUUID() :
	DbDataType()
{}

std::string
DbUUID::print(
	utility::sql_database::sessionOP db_session
) const {

	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		return "BLOB";
		break;
	case utility::sql_database::DatabaseMode::postgres :
		return "UUID";
		break;
	case utility::sql_database::DatabaseMode::mysql :
		//I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
		return "BINARY(16)";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

} // schema_generator
} // namespace database
} // namespace utility
