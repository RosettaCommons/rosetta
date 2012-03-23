// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BatchFeatures.cc
///
/// @brief
/// @author tim

// Unit Headers
#include <protocols/features/BatchFeatures.hh>

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

// C++ Headers
#include <string>
#include <sstream>

static basic::Tracer TR("protocols.features.BatchFeatures");

namespace protocols{
namespace features{

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

BatchFeatures::BatchFeatures( BatchFeatures const & ) :
FeaturesReporter()
{}

BatchFeatures::~BatchFeatures(){}

string
BatchFeatures::type_name() const { return "BatchFeatures"; }

string
BatchFeatures::schema() const {
	using namespace basic::database::schema_generator;


	PrimaryKey batch_id(
		Column("batch_id", DbInteger(), false /*not null*/, true /*autoincrement*/));

	ForeignKey protocol_id(
		Column("protocol_id", DbInteger()),
		"protocols",
		"protocol_id",
		true /*defer*/);

	Column name(Column("name", DbText()));

	Column description(Column("description", DbText()));


	Schema batches("batches", batch_id);
	batches.add_foreign_key(protocol_id);
	batches.add_column(name);
	batches.add_column(description);

	return batches.print();
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
	Size protocol_id,
	std::string name,
	std::string description,
	sessionOP db_session
){
	TR << "Writing to batches table with name " << name << ", referencing protocol " << protocol_id << std::endl;
	std::string insert_string("INSERT INTO batches (protocol_id, name, description) VALUES (?,?,?);");
	cppdb::statement insert_statement = basic::database::safely_prepare_statement(insert_string,db_session);
	insert_statement.bind(1,protocol_id);
	insert_statement.bind(2,name);
	insert_statement.bind(3,description);

	basic::database::safely_write_to_database(insert_statement);

	std::string select = "SELECT batch_id, protocol_id, name FROM batches;";
	cppdb::statement select_stmt = basic::database::safely_prepare_statement(select,db_session);

	cppdb::result res(basic::database::safely_read_from_database(select_stmt));
	while(res.next()){
		signed long long b_id, p_id;
		std::string test_name;
		res >> b_id >> p_id >> test_name;
	}

	int test = insert_statement.sequence_last("batches_batch_id_seq");
	return test;
}

} // features namesapce
} // protocols namespace
