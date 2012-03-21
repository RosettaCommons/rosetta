// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DbDataType.cc
///
/// @brief

/// @author tim

#include <utility/sql_database/DbDataType.hh>

//Basic
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Utility
#include <utility/exit.hh>
#include <utility/string_util.hh>

//C++
#include <string>

namespace utility{
namespace sql_database{

/*********Base class*********/
DbDataType::DbDataType() {
    if(basic::options::option[basic::options::OptionKeys::inout::database_mode].user()){
        database_mode_=basic::options::option[basic::options::OptionKeys::inout::database_mode].value();
    }
    else{
        database_mode_="sqlite3";
    }
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
    if(this->database_mode_.compare("sqlite3") == 0){
        type_string_ = "TEXT";
    }
    else if(this->database_mode_.compare("mysql") == 0 || this->database_mode_.compare("postgres") == 0){
        type_string_ = "VARCHAR(" + size_string + ")";
    }
    else{
        utility_exit_with_message("ERROR: Invalid database mode supplied. Please specify sqlite3, mysql, or postgres");
    }
    
}

//Needed because MYSQL doesn't support variable length primary keys
DbTextKey::DbTextKey():
DbDataType()
{
    
    if(this->database_mode_.compare("sqlite3") == 0){
        type_string_ = "TEXT";
    }
    else if(this->database_mode_.compare("postgres") == 0){
        type_string_ = "TEXT";
    }
    else if(this->database_mode_.compare("mysql") == 0){
        type_string_ = "VARCHAR(255)";
    }
    else{
        utility_exit_with_message("ERROR: Invalid database mode supplied. Please specify sqlite3, mysql, or postgres");
    }
    
}
    
/*********Integer data types*********/
DbInteger::DbInteger():
DbDataType()
{
    type_string_ = "INTEGER";
}

DbInteger::DbInteger(int size):
DbDataType()
{
    std::string size_string = utility::to_string(size);
    type_string_ = "INTEGER";
}
    
DbBigInt::DbBigInt():
DbDataType()
{    
    if(this->database_mode_.compare("sqlite3") == 0){
        type_string_ = "INTEGER";
    }
    else if(this->database_mode_.compare("mysql") == 0 || this->database_mode_.compare("postgres") == 0){
        type_string_ = "BIGINT";
    }
    else{
        utility_exit_with_message("ERROR: Invalid database mode supplied. Please specify sqlite3, mysql, or postgres");
    }
}
    
/*********Boolean data types*********/
DbBoolean::DbBoolean():
DbDataType()
{
    type_string_ = "BOOLEAN";
}
    
DbReal::DbReal():
DbDataType()
{
    type_string_ = "REAL";
}
    
/*********UUID data types*********/
DbUUID::DbUUID():
DbDataType()
{
    if(this->database_mode_.compare("sqlite3") == 0){
        type_string_ = "BLOB";
    }
    else if(this->database_mode_.compare("postgres") == 0){
        type_string_ = "UUID";
    }
    //I can't figure out how to write the 16byte UUID into a varbinary(16). I don't like mysql.
    else if(this->database_mode_.compare("mysql") == 0){
        type_string_ = "BINARY(36)";
    }
    else{
        utility_exit_with_message("ERROR: Invalid database mode supplied. Please specify sqlite3, mysql, or postgres");
    }
}
    
} // namespace sql_database
} // namespace utility
