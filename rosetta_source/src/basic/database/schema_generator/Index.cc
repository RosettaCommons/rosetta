// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Index.cc
///
/// @brief Index class for the schema generator framework
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <basic/database/schema_generator/Index.hh>
#include <basic/database/schema_generator/Column.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <platform/types.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/sql_database/types.hh>

//C++ Headers
#include <string>
#include <sstream>

static basic::Tracer TR("utility.sql_database.Index");

namespace basic{
namespace database{
namespace schema_generator{

using platform::Size;
using std::string;
using std::stringstream;
using utility::vector1;

Index::Index(
	Column column,
	bool unique
) :
	unique_(unique),
	database_mode_(),
	columns_()
{
	columns_.push_back(column);
	init_db_mode();
}

Index::Index(
	Columns columns,
	bool unique
) :
	unique_(unique),
	database_mode_(),
	columns_(columns)
{
	init_db_mode();
}

Index::Index(
	Index const & src
) :
	unique_(src.unique_),
	database_mode_(src.database_mode_),
	columns_(src.columns_)
{}

void
Index::init_db_mode(){
	database_mode_ =
		utility::sql_database::database_mode_from_name(
			basic::options::option[basic::options::OptionKeys::inout::dbms::mode]);
}

Columns
Index::columns(){
	return this->columns_;
}

string
Index::print(
	string const & table_name
){
	stringstream s;

	if(database_mode_ == utility::sql_database::DatabaseMode::sqlite3){
		s << "CREATE ";
		if(unique_){
			s << "UNIQUE ";
		}
		s << "INDEX IF NOT EXISTS\n\t";
		s << table_name;
		for(Size i=1; i <= columns_.size(); ++i){
			s << "_" << columns_[i].name();
		}
		s << " ON\n\t";
		s << " " << table_name << " ( ";
		for(Size i=1; i <= columns_.size(); ++i){
			if(i != 1){
				s << ", ";
			}
			s << columns_[i].name();
		}
		s << " );\n";

		return s.str();
	}

	return "";
}


} // schema_generator
} // namespace database
} // namespace utility
