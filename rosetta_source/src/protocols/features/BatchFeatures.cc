// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
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

static basic::Tracer TR("protocols.features.BatchFeatures");

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

void
BatchFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;

	PrimaryKey batch_id(
		Column("batch_id", new DbInteger(), false /*not null*/, false /*autoincrement*/));

	ForeignKey protocol_id(
		Column("protocol_id", new DbInteger()),
		"protocols",
		"protocol_id",
		true /*defer*/);

	Column name(Column("name", new DbText()));
	Column description(Column("description", new DbText()));

	Schema batches("batches", batch_id);
	batches.add_foreign_key(protocol_id);
	batches.add_column(name);
	batches.add_column(description);

	batches.write(db_session);
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
	TR.Debug
		<< "Writing to batches table with:" << std::endl
		<< "\tbatch_id '" << batch_id << "'" << std::endl
		<< "\tprotocol_id '" << protocol_id << "'" << std::endl
		<< "\tname '" << name << "'" << std::endl
		<< "\tdescription '" << description << "'" << std::endl;


	using namespace boost::assign;
	std::vector<string> column_names = list_of("batch_id")("protocol_id")("name")("description");

	stringstream batch_id_str;
	batch_id_str << batch_id;

	stringstream protocol_id_str;
	protocol_id_str << protocol_id;

	std::vector<string> values =
		list_of(batch_id_str.str())(protocol_id_str.str())("'"+name+"'")("'"+description+"'");
	basic::database::insert_or_ignore("batches", column_names, values, db_session);

	return 0;
}

} // features namesapce
} // protocols namespace
