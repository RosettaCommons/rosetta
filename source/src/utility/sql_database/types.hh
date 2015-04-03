// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/sql_database/types.hh
/// @author Matthew O'Meara


#ifndef INCLUDED_utility_sql_database_types_HH
#define INCLUDED_utility_sql_database_types_HH


// C++ Headers
#include <string>

namespace utility {
namespace sql_database {


struct TransactionMode {
	enum e {
		none = 1,
		standard,
		chunk
	};
};
	
struct DatabaseMode {
	enum e {
		sqlite3 = 1,
		mysql,
		postgres
	};
};

TransactionMode::e
transaction_mode_from_name(
  std::string transaction_mode);

std::string
name_from_transaction_mode(
  TransactionMode::e transaction_mode);
 
DatabaseMode::e
database_mode_from_name(
  std::string database_mode);

std::string
name_from_database_mode(
  DatabaseMode::e database_mode);


}
}


#endif
