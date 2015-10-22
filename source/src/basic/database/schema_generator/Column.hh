// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Column.hh
/// @brief Column class for the schema generator framework
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_schema_generator_Column_HH
#define INCLUDED_basic_database_schema_generator_Column_HH

// Unit Headers
#include <iosfwd>                                              // for string
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>  // for sessionOP
#include <basic/database/schema_generator/Column.fwd.hh>
#include <basic/database/schema_generator/DbDataType.fwd.hh>   // for DbData...
#include <platform/types.hh>                                   // for Size
#include <utility/pointer/ReferenceCount.hh>                   // for Refere...

namespace basic {
namespace database {
namespace schema_generator {

class Column : public utility::pointer::ReferenceCount {
public:

	Column(std::string name, DbDataTypeOP type);

	Column(std::string name, DbDataTypeOP type, bool allow_null);

	Column(std::string name, DbDataTypeOP type, bool allow_null, bool auto_increment, platform::Size auto_increment_base = 0);

	Column(Column const & src);

	virtual ~Column();

	std::string name() const;

	bool auto_increment() const;
	Size auto_increment_base() const;

	std::string print(utility::sql_database::sessionOP db_session) const;

	bool operator==(const Column &other) const;

private:

	std::string name_;
	DbDataTypeOP type_;
	bool allow_null_;
	bool auto_increment_;
	platform::Size auto_increment_base_;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
