// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Column.cc
///
/// @brief Column class for the schema generator framework
/// @author Tim Jacobs

//Unit Headers
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/sql_database/types.hh>

namespace basic{
namespace database{
namespace schema_generator{

Column::Column(std::string name, DbDataType type) :
	name_(name),
	type_(type),
	allow_null_(true),
	auto_increment_(false)
{
	init_db_mode();
}

Column::Column(std::string name, DbDataType type, bool allow_null) :
	name_(name),
	type_(type),
	allow_null_(allow_null),
	auto_increment_(false)
{
	init_db_mode();
}

Column::Column(std::string name, DbDataType type, bool allow_null, bool auto_increment) :
	name_(name),
	type_(type),
	allow_null_(allow_null),
	auto_increment_(auto_increment)
{
	init_db_mode();
}

Column::Column(Column const & src) :
	ReferenceCount(),
	name_(src.name_),
	type_(src.type_),
	allow_null_(src.allow_null_),
	auto_increment_(src.auto_increment_)
{
	init_db_mode();
}

Column::~Column() {}

void Column::init_db_mode(){
	database_mode_ =
		utility::sql_database::database_mode_from_name(
			basic::options::option[basic::options::OptionKeys::inout::dbms::mode]);
}

std::string Column::name() const{
	return name_;
}

bool Column::auto_increment() const{
	return this->auto_increment_;
}

std::string Column::print() const{
	std::string column_string = "";
	if(auto_increment_){
		column_string += name_ + " ";
		switch(database_mode_) {
		case utility::sql_database::DatabaseMode::sqlite3:
			column_string += this->type_.print() + " PRIMARY KEY AUTOINCREMENT"; //only way to autoincrement in SQLite is with a primary key
			name_ + " " + type_.print();
			break;
		case utility::sql_database::DatabaseMode::mysql:
			column_string += this->type_.print() + " AUTO_INCREMENT";
			break;
		case utility::sql_database::DatabaseMode::postgres:
			column_string += "BIGSERIAL";
			break;
		default:
			utility_exit_with_message("ERROR:Please specify the database mode using -inout::dbms::mode. Valid options are: 'sqlite3', 'mysql', or 'postgres'");
		}
	}
	else{
		column_string += this->name_ + " " + this->type_.print();
	}

	if(!allow_null_){
		column_string += " NOT NULL";
	}
	return column_string;
}

bool Column::operator==(const Column &other) const {
	return (this->name_.compare(other.name()) == 0);
}

} // schema_generator
} // namespace database
} // namespace utility
