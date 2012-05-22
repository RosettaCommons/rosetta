// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PrimaryKey.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_basic_database_schema_generator_PrimaryKey_HH
#define INCLUDED_basic_database_schema_generator_PrimaryKey_HH

#include <basic/database/schema_generator/Column.hh>
#include <utility/vector1.hh>

//C++ Header
#include <string>
#include <set>

namespace basic{
namespace database{
namespace schema_generator{

class PrimaryKey
{
public:

	PrimaryKey();

	PrimaryKey(Column column);

	PrimaryKey(utility::vector1<Column> columns);

	void add_column(Column column);

	utility::vector1<Column> columns();

	std::string print();

private:
	utility::vector1<Column> columns_;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
