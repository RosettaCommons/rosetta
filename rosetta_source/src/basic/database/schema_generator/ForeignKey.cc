// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ForeignKey.cc
///
/// @brief
/// @author tim

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

ForeignKey::ForeignKey(Column column, std::string reference_table, std::string reference_column):
column_(column),
reference_table_(reference_table),
reference_column_(reference_column),
defer_(false)
{
	init_db_mode();
}

ForeignKey::ForeignKey(Column column, std::string reference_table, std::string reference_column, bool defer):
column_(column),
reference_table_(reference_table),
reference_column_(reference_column),
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

std::string ForeignKey::print(){
	std::string foreign_key_string = "FOREIGN KEY (" + column_.name() + ") REFERENCES " + reference_table_ + "(" + reference_column_ + ")";
	if(defer_){

		if(this->database_mode_.compare("sqlite3") == 0 || this->database_mode_.compare("postgres") == 0){
			foreign_key_string += " DEFERRABLE INITIALLY DEFERRED";
		}
		else if(this->database_mode_.compare("mysql") == 0){
			//MySQL does not support deferring foreign keys. Warn and continue
			TR << "Warning: You have tried to defer and foreign key constraint in MySql mode. MySql does not support foreign keys!" << std::endl;
		}
		else{
			utility_exit_with_message("ERROR:Please specify the database mode using -inout::database_mode. Valid options are: 'sqlite3', 'mysql', or 'postgres'");
		}
	}
	return foreign_key_string;
}

Column ForeignKey::column(){
	return this->column_;
}

} // schema_generator
} // namespace database
} // namespace utility

