// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/sql_database/DatabaseSessionManagerTests.cxxtest.hh
/// @brief  Test DatabaseSessionManager
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>

// Unit Headers
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <iostream>

// Can't use tracers because we're in utility
using std::cout;
using std::cerr;
using std::endl;

class DatabaseSessionManagerTests : public CxxTest::TestSuite {

public:

  void
  setup() {}

  void test_cppdb_interface() {

    cppdb::session db_session("sqlite3:db=test_cppdb_interface.db3");
    TS_ASSERT(true);

  }

  void test_db_session(){

    using namespace utility::sql_database;
    using namespace cppdb;
    using std::string;

    DatabaseSessionManager * scm( DatabaseSessionManager::get_instance() );

    try {
    	sessionOP db_session( scm->get_session_sqlite3("test_db_session.db3"));

      cppdb::statement drop_table = (*db_session) << "DROP TABLE IF EXISTS users";
    	drop_table.exec();
      cppdb::statement create_table = (*db_session) <<
        "CREATE TABLE users ( "
        " id integer primary key not null, "
        " name varchar(128) not null);";
    	create_table.exec();
      cppdb::statement insert_Moshe = (*db_session) <<
        "INSERT INTO users(id,name) VALUES(?,?)" << 1 << "Moshe";
      insert_Moshe.exec();
      cppdb::statement insert_Yossi = (*db_session) <<
        "INSERT INTO users(id,name) VALUES(?,?)" << 2 << "Yossi";
      insert_Yossi.exec();
      cppdb::result name1 = (*db_session) <<
        "SELECT name FROM users WHERE id=?" << 1;
      if( name1.next() ){
        string name;
        name1 >> name;
        //cout<<name<<endl;
    	} else {
    	  //cout<<"No user with id="<<1<<endl;
    	}
      cppdb::result id_names = (*db_session) << "SELECT id,name FROM users";
    	while(id_names.next()) {
    	  int id;
    	  string name;
    	  id_names >> id >> name;
    	  //cout << id << "\t" << name << endl;
    	}
    } catch(std::exception const &e) {
      cerr << e.what() << endl;
      TS_ASSERT(false);
    }
  }

  void test_multiple_selections() {

    using namespace utility::sql_database;
    using namespace cppdb;
    using std::string;

    DatabaseSessionManager * scm( DatabaseSessionManager::get_instance() );

    sessionOP db_session1(scm->get_session_sqlite3("test_db_session.db3"));
    sessionOP db_session2(scm->get_session_sqlite3("test_cppdb_interface.db3"));


  }


};


