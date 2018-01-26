// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/DatabaseResourceLocatorTests.cxxtest.hh
/// @brief  Test the DatabaseResourceLocator class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>

#include <basic/Tracer.hh>

#include <basic/resource_manager/locator/DatabaseResourceLocator.hh>
#include <basic/resource_manager/locator/StringResourceStream.hh>

#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ headers
#include <cstdio>
#include <string>
#include <sstream>

static basic::Tracer tr( "basic.resource_manager.locator.DatabaseResourceLocator.cxxtest");

class DatabaseResourceLocatorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_DatabaseResourceLocator_query_sqlite_database() {
		using namespace basic::resource_manager;
		using namespace basic::resource_manager::locator;
		using namespace utility::sql_database;
		using namespace cppdb;

		std::string test_db_name( "DatabaseResourceLocator.db3" );

		{
			// Create the Database that we'll be reading from later.

			// If we ran this unit test before, then the DB might already be there: delete it.
			if ( FILE * file = fopen( test_db_name.c_str(), "r" ) ) {
				fclose( file );
				remove( test_db_name.c_str() );
			}

			sessionOP db_session = basic::database::get_db_session( test_db_name );
			using namespace basic::database::schema_generator;

			Column id( "id", DbDataTypeOP( new DbBigInt ) );
			Column value( "val", DbDataTypeOP( new DbText ) );
			PrimaryKey primary_key( id );
			Schema table( "example_table", primary_key );
			table.add_column( value );
			table.write( db_session );

			std::string statement_string = "INSERT INTO example_table (id, val) VALUES (?,?)";
			statement stmt(basic::database::safely_prepare_statement(statement_string, db_session ));
			std::list< std::pair< int, std::string > > vals_to_insert = {
				{ 1, "Unit"},
				{ 2, "Test"},
				{ 3, "Your"},
				{ 4, "Code"}
				};

			for ( auto const & val_pair : vals_to_insert ) {
				stmt.bind(1,val_pair.first);
				stmt.bind(2,val_pair.second);
				basic::database::safely_write_to_database( stmt );
			}
		}

		std::string db_locator_tag_string =
			"<DatabaseResourceLocator database_name=\"" +
			test_db_name + "\" sql_command=\"SELECT val FROM example_table WHERE id = ?;\"/>";

		utility::tag::TagCOP db_locator_tag = utility::tag::Tag::create( db_locator_tag_string );;
		DatabaseResourceLocator db_locator;
		db_locator.parse_my_tag( db_locator_tag );

		ResourceStreamOP rstream = db_locator.locate_resource_stream( "3" );
		std::string test_contents;
		std::getline( rstream->stream(), test_contents );
		TS_ASSERT_EQUALS( test_contents, "Your" );

		TS_ASSERT( dynamic_cast< StringResourceStream * > ( rstream.get() ) );
	}


};
