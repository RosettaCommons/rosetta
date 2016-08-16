// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/database/schema_generator/PrimaryKey.hh
///
/// @brief PrimaryKey class for the schema generator framework
/// @author Tim Jacobs


#ifndef INCLUDED_basic_database_schema_generator_PrimaryKey_HH
#define INCLUDED_basic_database_schema_generator_PrimaryKey_HH

#include <basic/database/schema_generator/Column.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

//C++ Header
#include <iosfwd>

namespace basic {
namespace database {
namespace schema_generator {

class PrimaryKey
{
public:

	PrimaryKey();

	PrimaryKey(Column column);

	PrimaryKey(Columns columns);

	void add_column(Column column);

	Columns const & columns() const;

	std::string
	print(
		utility::sql_database::sessionOP db_session
	) const;

private:
	Columns columns_;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
