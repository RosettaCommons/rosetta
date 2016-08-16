// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Schema.hh
///
/// @brief table definition for the schema generator framework
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_schema_generator_Schema_HH
#define INCLUDED_basic_database_schema_generator_Schema_HH

#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Index.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <iosfwd>
#include <string>

#include <basic/database/schema_generator/Column.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.hh>


namespace basic {
namespace database {
namespace schema_generator {

class Schema
{

public:

	Schema(
		std::string table_name);

	Schema(
		std::string table_name,
		PrimaryKey primary_key);

	Schema(
		Schema const & src);

	void add_foreign_key(ForeignKey key);

	void add_column(Column column);

	void add_constraint(ConstraintOP constraint);

	void add_index(Index index);

	//@brief Returns sql statements needed to fully initialize the given schema.
	//
	//See table_schema_statements and table_init_statements, statement blocks
	//must be executed separately in sqlite.
	std::string print(utility::sql_database::sessionOP db_session) const;

	//@brief Write schema to database and initialize if needed.
	//
	//Retries write once to resolve transaction conflicts.
	void write(utility::sql_database::sessionOP db_session);

protected:
	//@brief Returns sql statements to declare table, keys, and indices.
	std::string table_schema_statements(utility::sql_database::sessionOP db_session) const;

	//@brief Returns sql to initialize table state.
	//
	//Currently limited to autoincrement starting values. Must be executed in separate statement
	//when using sqlite3, as sqlite_sequence table isn't instantiated until the statement containing
	//CREATE TABLE with autoincrement value is executed.
	std::string table_init_statements(utility::sql_database::sessionOP db_session) const;

	//@brief Write and initialize schema if not already declared.
	void check_table_and_perform_write(
		utility::sql_database::sessionOP db_session,
		std::string schema_statement,
		std::string init_statements) const;

private:

	void init();

	std::string table_name_;
	PrimaryKey primary_key_;
	Columns columns_;
	utility::vector1<ForeignKey> foreign_keys_;
	utility::vector1<ConstraintOP> constraints_;
	utility::vector1<Index> indices_;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
