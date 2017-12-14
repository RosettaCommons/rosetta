// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/schema_generator/ForeignKey.cc
///
/// @brief ForeignKey class for the schema generator framework
/// @author Tim Jacobs

#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <utility>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic Headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "utility.sql_database.ForeignKey" );


// Utility Headers
#include <utility/exit.hh>

//C++ Headers
#include <string>

namespace basic {
namespace database {
namespace schema_generator {

using std::string;
using utility::vector1;

ForeignKey::ForeignKey(
	Column const & column,
	std::string const & reference_table,
	std::string const & reference_column) :
	columns_(),
	reference_columns_(),
	reference_table_(reference_table),
	defer_(false)
{
	columns_.push_back(column);
	reference_columns_.push_back(reference_column);
}

ForeignKey::ForeignKey(
	Column const & column,
	string const & reference_table,
	string const & reference_column,
	bool defer) :
	columns_(),
	reference_columns_(),
	reference_table_(reference_table),
	defer_(defer)
{
	columns_.push_back(column);
	reference_columns_.push_back(reference_column);
}

ForeignKey::ForeignKey(
	Columns const & columns,
	string const & reference_table,
	vector1<string> const & reference_columns,
	bool defer) :
	columns_(columns),
	reference_columns_(reference_columns),
	reference_table_(reference_table),
	defer_(defer)
{}

Columns ForeignKey::columns(){
	return this->columns_;
}

std::string
ForeignKey::print(
	utility::sql_database::sessionOP db_session
) const {
	std::string foreign_key_string = "FOREIGN KEY (";

	for ( size_t i=1; i<=columns_.size(); ++i ) {
		foreign_key_string += columns_[i].name();
		if ( i != columns_.size() ) {
			foreign_key_string+=", ";
		}
	}
	foreign_key_string += ") REFERENCES " + reference_table_ + "(";

	for ( size_t i=1; i<=reference_columns_.size(); ++i ) {
		foreign_key_string += reference_columns_[i];
		if ( i != reference_columns_.size() ) {
			foreign_key_string+=", ";
		}
	}
	foreign_key_string += ")";

	if ( defer_ ) {
		switch(db_session->get_db_mode()) {
		case utility::sql_database::DatabaseMode::mysql :
			//MySQL does not support deferring foreign keys.
			break;
		case utility::sql_database::DatabaseMode::postgres :
			foreign_key_string += " DEFERRABLE INITIALLY DEFERRED";
			break;
		case utility::sql_database::DatabaseMode::sqlite3 :
			foreign_key_string += " DEFERRABLE INITIALLY DEFERRED";
			break;
		default :
			utility_exit_with_message(
				"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
		}
	}
	return foreign_key_string;
}

} // schema_generator
} // namespace database
} // namespace utility
