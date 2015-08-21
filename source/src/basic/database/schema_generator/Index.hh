// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Index.hh
/// @brief Index class for the schema generator framework
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_basic_database_schema_generator_Index_HH
#define INCLUDED_basic_database_schema_generator_Index_HH

#include <utility/vector1.hh>
#include <utility/sql_database/types.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <basic/database/schema_generator/Column.hh>

//C++ Headers
#include <string>

namespace basic {
namespace database {
namespace schema_generator {

class Index {
public:

	Index(
		Column column,
		bool unique=true);

	Index(
		Columns columns,
		bool unique=true);

	Index(
		Index const & src);

	Columns
	columns();

	bool
	unique() { return unique_; }

	std::string
	print(
		std::string const & table_name,
		utility::sql_database::sessionOP db_session
	) const ;

private:
	bool unique_;
	Columns columns_;
};


} // schema_generator
} // namespace database
} // namespace utility

#endif
