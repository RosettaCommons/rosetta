// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/schema_generator/Index.cc
///
/// @brief Index class for the schema generator framework
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <basic/database/schema_generator/Index.hh>
#include <basic/database/schema_generator/Column.hh>
#include <utility>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <platform/types.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/sql_database/types.hh>

//C++ Headers
#include <string>
#include <sstream>

static basic::Tracer TR( "utility.sql_database.Index" );

namespace basic {
namespace database {
namespace schema_generator {

using platform::Size;
using std::string;
using std::stringstream;
using utility::vector1;

Index::Index(
	Column column,
	bool unique
) :
	unique_(unique),
	columns_()
{
	columns_.push_back(column);
}

Index::Index(
	Columns columns,
	bool unique
) :
	unique_(unique),
	columns_(std::move(columns))
{}

Index::Index( Index const & /*src*/ ) = default;

Columns
Index::columns(){
	return this->columns_;
}

string
Index::print(
	string const & table_name,
	utility::sql_database::sessionOP db_session
) const {
	stringstream s;

	switch(db_session->get_db_mode()) {
	case utility::sql_database::DatabaseMode::mysql :
		break;
	case utility::sql_database::DatabaseMode::postgres :
		break;
	case utility::sql_database::DatabaseMode::sqlite3 : {
		s << "CREATE ";
		if ( unique_ ) {
			s << "UNIQUE ";
		}
		s << "INDEX IF NOT EXISTS\n\t";
		s << table_name;
		for ( Size i=1; i <= columns_.size(); ++i ) {
			s << "_" << columns_[i].name();
		}
		s << " ON\n\t";
		s << " " << table_name << " ( ";
		for ( Size i=1; i <= columns_.size(); ++i ) {
			if ( i != 1 ) {
				s << ", ";
			}
			s << columns_[i].name();
		}
		s << " );\n";

		break;
	}
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_session->get_db_mode()) + "'");
		return "";  // just here to remove warning; lame ~Labonte
	}

	return s.str();
}


} // schema_generator
} // namespace database
} // namespace utility
