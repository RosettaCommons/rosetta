// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   utility/sql_database/types.hh
/// @author Matthew O'Meara


#ifndef INCLUDED_utility_sql_database_types_HH
#define INCLUDED_utility_sql_database_types_HH


// C++ Headers
#include <string>

namespace utility {
namespace sql_database {



struct DatabaseMode {
	enum e {
		sqlite3 = 1,
		mysql,
		postgres
	};
};

DatabaseMode::e
database_mode_from_name(
  std::string database_mode);

std::string
name_from_database_mode(
  DatabaseMode::e database_mode);



}
}


#endif
