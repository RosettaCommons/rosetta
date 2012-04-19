// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Schema.cc
///
/// @brief Construct a database backend independant schema
/// @author Tim Jacobs

//Unit
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility Headers
#include <utility/exit.hh>

#include <string>
#include <stdio.h>
#include <set>

namespace basic{
namespace database{
namespace schema_generator{

Schema::Schema(std::string table_name):
table_name_(table_name)
{
	init();
}

Schema::Schema(std::string table_name, PrimaryKey primary_key):
table_name_(table_name),
primary_key_(primary_key)
{
	init();
}

void Schema::init(){
	if(basic::options::option[basic::options::OptionKeys::inout::database_mode].user()){
		database_mode_=basic::options::option[basic::options::OptionKeys::inout::database_mode].value();
	}
	else{
		database_mode_="sqlite3";
	}

	//Add primary key columns to schema list
	utility::vector1<Column> key_columns = primary_key_.columns();
	this->columns_.insert( columns_.end(), key_columns.begin(), key_columns.end() );

}

void Schema::add_foreign_key(ForeignKey key){
	this->foreign_keys_.push_back(key);
	//if the foreign key is also a primary key it will have already been added

	utility::vector1<Column> key_cols = key.columns();

	for(size_t i=1; i <= key_cols.size(); ++i){
		if(!this->columns_.contains(key_cols[i]))
		{
			this->columns_.push_back(key_cols[i]);
		}
	}
}

void Schema::add_column(Column column){
	//Don't add a column more than once
	if(!this->columns_.contains(column))
	{
		this->columns_.push_back(column);
	}
}

void Schema::add_constraint(ConstraintOP constraint){
	this->constraints_.push_back(constraint);
}

std::string Schema::print(){
	std::string schema_string = "CREATE TABLE IF NOT EXISTS " + table_name_ + "(\n\t";

	for (utility::vector1<Column>::const_iterator it=columns_.begin(); it!=columns_.end(); it++){
		if(it!=columns_.begin()){
			schema_string += ",\n\t";
		}
		schema_string += it->print();
	}

	for(size_t i=1; i<=foreign_keys_.size(); i++){
		schema_string += ",\n\t" + foreign_keys_[i].print();
	}

	if(primary_key_.columns().size() > 0){
		if(database_mode_ != "sqlite3"){
			schema_string += ",\n\t" + primary_key_.print();
		}
		else{
			//Prevent adding the primary key twice - this will happen if you have an autoincrementing primary key in sqlite3
			utility::vector1<Column> keys = this->primary_key_.columns();

			if(!(keys.size()==1 && keys.begin()->auto_increment())){
				schema_string += ",\n\t" + primary_key_.print();
			}
		}
	}

	for(size_t i=1; i<=constraints_.size(); i++){
		schema_string += ",\n\t" + constraints_[i]->print();
	}

	schema_string += ");";
	return schema_string;
}

} // schema_generator
} // namespace database
} // namespace utility
