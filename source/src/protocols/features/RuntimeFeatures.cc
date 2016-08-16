// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
// Unit headers
#include <protocols/features/RuntimeFeatures.hh>
#include <protocols/features/RuntimeFeaturesCreator.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// External headers
#include <cppdb/frontend.h>

// C++ headers
#include <string>

// Namespaces {{{1
using namespace std;
using core::Size;
using core::Real;
using core::pose::Pose;
using utility::vector1;
using utility::tools::make_vector1;
using utility::sql_database::sessionOP;
// }}}1

namespace protocols {
namespace features {

FeaturesReporterOP RuntimeFeaturesCreator::create_features_reporter() const { // {{{1
	return FeaturesReporterOP( new RuntimeFeatures );
}

string RuntimeFeaturesCreator::type_name() const { // {{{1
	return "RuntimeFeatures";
}
// }}}1

RuntimeFeatures::RuntimeFeatures() {} // {{{1

RuntimeFeatures::~RuntimeFeatures() {} // {{{1

vector1 <string> RuntimeFeatures::features_reporter_dependencies() const { // {{{1
	return make_vector1("StructureFeatures");
}

void RuntimeFeatures::write_schema_to_db(sessionOP db_session) const { // {{{1
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt ));
	Column timestamp("timestamp", DbDataTypeOP( new DbText(20) ));
	Column elapsed_time("elapsed_time", DbDataTypeOP( new DbInteger ));

	PrimaryKey primary_key(struct_id);
	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);

	Schema table("runtimes", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(timestamp);
	table.add_column(elapsed_time);
	table.write(db_session);
}

Size RuntimeFeatures::report_features( // {{{1
	Pose const & /*pose*/,
	vector1 <bool> const & /*relevant_residues*/,
	StructureID struct_id,
	sessionOP db_session) {

	using protocols::jd2::JobCOP;
	using protocols::jd2::JobDistributor;

	JobCOP job = JobDistributor::get_instance()->current_job();

	string statement_string =
		"INSERT INTO runtimes (struct_id, timestamp, elapsed_time) VALUES (?,?,?);";
	cppdb::statement statement = basic::database::safely_prepare_statement(
		statement_string, db_session);

	statement.bind(1, struct_id);
	statement.bind(2, job->timestamp());
	statement.bind(3, job->elapsed_time());

	basic::database::safely_write_to_database(statement);
	return 0;
}
// }}}1

}
}
