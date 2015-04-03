// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/BatchFeatures.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/features/BatchFeatures.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>

// Project Headers
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <core/types.hh>
#include <core/svn_version.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

//Basic Headers
#include <basic/Tracer.hh>

// External Headers
#include <cppdb/frontend.h>
#include <boost/assign/list_of.hpp>

// C++ Headers
#include <string>
#include <sstream>


namespace protocols{
namespace features{

static thread_local basic::Tracer TR( "protocols.features.BatchFeatures" );

using std::string;
using std::stringstream;
using basic::options::OptionKeys::parser::protocol;
using basic::options::option;
using core::Size;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

BatchFeatures::BatchFeatures(){}

BatchFeatures::BatchFeatures( BatchFeatures const & ) : utility::pointer::ReferenceCount() {}

BatchFeatures::~BatchFeatures(){}

string
BatchFeatures::type_name() const { return "BatchFeatures"; }

void
BatchFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session, core::Size batch_id) const{

	using namespace basic::database::schema_generator;

	if(batch_id){
		PrimaryKey batch_id(Column("batch_id", DbDataTypeOP( new DbInteger() )));

		ForeignKey protocol_id(
			Column("protocol_id", DbDataTypeOP( new DbInteger() )),
			"protocols",
			"protocol_id",
			true /*defer*/);

		Column name(Column("name", DbDataTypeOP( new DbText() )));
		Column description(Column("description", DbDataTypeOP( new DbText() )));

		Schema batches("batches", batch_id);
		batches.add_foreign_key(protocol_id);
		batches.add_column(name);
		batches.add_column(description);

		batches.write(db_session);

	}
	else{
		PrimaryKey batch_id(
			Column("batch_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, true /*autoincrement*/));

		ForeignKey protocol_id(
			Column("protocol_id", DbDataTypeOP( new DbInteger() )),
			"protocols",
			"protocol_id",
			true /*defer*/);

		Column name(Column("name", DbDataTypeOP( new DbText() )));
		Column description(Column("description", DbDataTypeOP( new DbText() )));

		Schema batches("batches", batch_id);
		batches.add_foreign_key(protocol_id);
		batches.add_column(name);
		batches.add_column(description);

		batches.write(db_session);
	}
}

utility::vector1<std::string>
BatchFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	return dependencies;
}


string
BatchFeatures::indices() const {
	return "";
}

Size
BatchFeatures::report_features(
	Size batch_id,
	Size protocol_id,
	std::string name,
	std::string description,
	sessionOP db_session
){
	cppdb::statement insert_statement;
	if(batch_id){

		//Check to make sure an existing batch with the same id doesn't exist
		std::string statement_string =
			"SELECT\n"
			"	count(*)\n"
			"FROM\n"
			"	batches\n"
			"WHERE\n"
			"	batch_id = ?;";
		cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,protocol_id);

		TR << "Checking for existing batches entry with given id" << std::endl;
		cppdb::result res(basic::database::safely_read_from_database(stmt));
		if(res.next()) {
			core::Size selected = 0;
			res >> selected;
			if(selected != 0) {
				return protocol_id;
			}
		}

		TR << "Writing to batches table with given batch id: " << batch_id << std::endl;
		std::string insert_string("INSERT INTO batches (batch_id, protocol_id, name, description) VALUES (?,?,?,?);");
		insert_statement = basic::database::safely_prepare_statement(insert_string,db_session);
		insert_statement.bind(1,batch_id);
		insert_statement.bind(2,protocol_id);
		insert_statement.bind(3,name);
		insert_statement.bind(4,description);

	} else {
		TR << "No batch ID, generating one automagically" << std::endl;
		std::string insert_string("INSERT INTO batches (protocol_id, name, description) VALUES (?,?,?);");
		insert_statement = basic::database::safely_prepare_statement(insert_string,db_session);
		insert_statement.bind(1,protocol_id);
		insert_statement.bind(2,name);
		insert_statement.bind(3,description);
	}

	//Batch features doesn't use safely_write_to_database due to special handling of thrown cppdb exceptions
	//basic::database::safely_write_to_database(insert_statement);
	insert_statement.exec();
	if(batch_id) {
		return batch_id;
	} else {
		core::Size new_batch_id = insert_statement.sequence_last("batches_batch_id_seq");
		return new_batch_id;
	}

//	TR.Debug
//		<< "Writing to batches table with:" << std::endl
//		<< "\tbatch_id '" << batch_id << "'" << std::endl
//		<< "\tprotocol_id '" << protocol_id << "'" << std::endl
//		<< "\tname '" << name << "'" << std::endl
//		<< "\tdescription '" << description << "'" << std::endl;

	return 0;
}

} // features namesapce
} // protocols namespace
