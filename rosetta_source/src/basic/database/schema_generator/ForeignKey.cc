// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/ForeignKey.cc
///
/// @brief ForeignKey class for the schema generator framework
/// @author Tim Jacobs

#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

static basic::Tracer TR("utility.sql_database.ForeignKey");


// Utility Headers
#include <utility/exit.hh>

//C++ Headers
#include <string>

namespace basic{
namespace database{
namespace schema_generator{

using std::string;
using utility::vector1;

ForeignKey::ForeignKey(
	Column column,
	std::string reference_table,
	std::string reference_column) :
	database_mode_(),
	columns_(),
	reference_columns_(),
	reference_table_(reference_table),
	defer_(false)
{
	columns_.push_back(column);
	reference_columns_.push_back(reference_column);
	init_db_mode();
}

ForeignKey::ForeignKey(
	Column column,
	string reference_table,
	string reference_column,
	bool defer) :
	database_mode_(),
	columns_(),
	reference_columns_(),
	reference_table_(reference_table),
	defer_(defer)
{
	columns_.push_back(column);
	reference_columns_.push_back(reference_column);
	init_db_mode();
}

ForeignKey::ForeignKey(
	Columns columns,
	string reference_table,
	vector1<string> reference_columns,
	bool defer) :
	database_mode_(),
	columns_(columns),
	reference_columns_(reference_columns),
	reference_table_(reference_table),
	defer_(defer)
{
	init_db_mode();
}

void ForeignKey::init_db_mode(){
	if(basic::options::option[basic::options::OptionKeys::inout::database_mode].user()){
		database_mode_=basic::options::option[basic::options::OptionKeys::inout::database_mode].value();
	}
	else{
		database_mode_="sqlite3";
	}
}

Columns ForeignKey::columns(){
	return this->columns_;
}

std::string ForeignKey::print(){
	std::string foreign_key_string = "FOREIGN KEY (";

	for(size_t i=1; i<=columns_.size(); ++i){
		foreign_key_string += columns_[i].name();
		if(i != columns_.size()){
			foreign_key_string+=", ";
		}
	}
	foreign_key_string += ") REFERENCES " + reference_table_ + "(";

	for(size_t i=1; i<=reference_columns_.size(); ++i){
		foreign_key_string += reference_columns_[i];
		if(i != reference_columns_.size()){
			foreign_key_string+=", ";
		}
	}
	foreign_key_string += ")";

	if(defer_){

		if(this->database_mode_.compare("sqlite3") == 0 || this->database_mode_.compare("postgres") == 0){
			foreign_key_string += " DEFERRABLE INITIALLY DEFERRED";
		}
		else if(this->database_mode_.compare("mysql") == 0){
			//MySQL does not support deferring foreign keys.
		}
		else{
			utility_exit_with_message("ERROR:Please specify the database mode using -inout::database_mode. Valid options are: 'sqlite3', 'mysql', or 'postgres'");
		}
	}
	return foreign_key_string;
}

} // schema_generator
} // namespace database
} // namespace utility
