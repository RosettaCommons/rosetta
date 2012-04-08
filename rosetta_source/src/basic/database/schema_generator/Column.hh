// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Column.hh
///
/// @brief
/// @author tim

#ifndef INCLUDED_basic_database_schema_generator_Column_HH
#define INCLUDED_basic_database_schema_generator_Column_HH

#include <basic/database/schema_generator/DbDataType.hh>

//C++ headers
#include <string>

namespace basic{
namespace database{
namespace schema_generator{

class Column
{
public:

	Column(std::string name, DbDataType type);

	Column(std::string name, DbDataType type, bool allow_null);
	
	Column(std::string name, DbDataType type, bool allow_null, bool auto_increment);

	void init_db_mode();

	std::string name() const;

	bool auto_increment() const;

	std::string print() const;

	bool operator==(const Column &other) const;

private:

	std::string database_mode_;
	std::string name_;
	DbDataType type_;
	bool allow_null_;
	bool auto_increment_;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
