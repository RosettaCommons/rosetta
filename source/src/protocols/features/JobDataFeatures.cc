// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/JobDataFeatures.cc
/// @author Sam DeLuca

//unit headers
#include <protocols/features/JobDataFeatures.hh>


//platform headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>


//External

//external headers
#include <cppdb/frontend.h>

#include <string>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>


namespace protocols {
namespace features {

using basic::database::table_exists;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

JobDataFeatures::JobDataFeatures() {}

JobDataFeatures::JobDataFeatures(JobDataFeatures const & ) : protocols::features::FeaturesReporter() {}

JobDataFeatures::~JobDataFeatures() {}

std::string JobDataFeatures::type_name() const
{
	return "JobDataFeatures";
}

void
JobDataFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false /*not null*/, false /*don't autoincrement*/);
	Column data_key("data_key", DbDataTypeOP( new DbText(255) ));

	utility::vector1<Column> primary_columns;
	primary_columns.push_back(struct_id);
	primary_columns.push_back(data_key);
	PrimaryKey primary_key(primary_columns);

	Schema job_string_data("job_string_data",primary_key);
	job_string_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	job_string_data.write(db_session);

	Schema job_string_string_data("job_string_string_data",primary_key);
	job_string_string_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	job_string_string_data.add_column(Column("data_value", DbDataTypeOP( new DbText() )));

	job_string_string_data.write(db_session);

	Schema job_string_real_data("job_string_real_data",primary_key);
	job_string_real_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	job_string_real_data.add_column(Column("data_value", DbDataTypeOP( new DbReal() )));

	job_string_real_data.write(db_session);

}

utility::vector1<std::string>
JobDataFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}


core::Size
JobDataFeatures::report_features(
	core::pose::Pose const & /*pose */,
	utility::vector1<bool> const & /*relevant_residues*/,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
)
{
	protocols::jd2::JobCOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	insert_string_rows(struct_id,db_session,job);
	insert_string_string_rows(struct_id,db_session,job);
	insert_string_real_rows(struct_id,db_session,job);
	return 0;
}

void
JobDataFeatures::load_into_pose(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::pose::Pose & pose
){
	load_string_data(db_session, struct_id, pose);
	load_string_string_data(db_session, struct_id, pose);
	load_string_real_data(db_session, struct_id, pose);
}

void JobDataFeatures::delete_record(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
)
{

	std::string delete_js_string = "DELETE FROM job_string_data WHERE struct_id = ?;\n";
	cppdb::statement delete_js_statement(basic::database::safely_prepare_statement(delete_js_string,db_session));
	delete_js_statement.bind(1,struct_id);
	basic::database::safely_write_to_database(delete_js_statement);

	std::string delete_ss_string = "DELETE FROM job_string_string_data WHERE struct_id = ?;\n";
	cppdb::statement delete_ss_statement(basic::database::safely_prepare_statement(delete_ss_string,db_session));
	delete_ss_statement.bind(1,struct_id);
	basic::database::safely_write_to_database(delete_ss_statement);

	std::string delete_sr_string = "DELETE FROM job_string_real_data WHERE struct_id = ?;";
	cppdb::statement delete_sr_statement(basic::database::safely_prepare_statement(delete_sr_string,db_session));
	delete_sr_statement.bind(1,struct_id);
	basic::database::safely_write_to_database(delete_sr_statement);

}

void JobDataFeatures::insert_string_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{

	InsertGenerator string_insert("job_string_data");
	string_insert.add_column("struct_id");
	string_insert.add_column("data_key");

	protocols::jd2::Job::Strings::const_iterator it(job->output_strings_begin());

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );
	for ( ; it != job->output_strings_end(); ++it ) {

		RowDataBaseOP string_data( new RowData<std::string>("data_key",*it) );

		string_insert.add_row(utility::tools::make_vector(struct_id_data,string_data));
	}

	string_insert.write_to_database(db_session);

}

void
JobDataFeatures::load_string_data(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::pose::Pose &
){
	if ( !table_exists(db_session, "job_string_data") ) return;

	std::string statement_string =
		"SELECT\n"
		"\tdata_key\n"
		"FROM\n"
		"\tjob_string_data\n"
		"WHERE\n"
		"\tjob_string_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while ( res.next() ) {
		std::string data_key;
		res >> data_key;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string(data_key);
	}
}

void JobDataFeatures::insert_string_string_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{

	InsertGenerator string_string_insert("job_string_string_data");
	string_string_insert.add_column("struct_id");
	string_string_insert.add_column("data_key");
	string_string_insert.add_column("data_value");

	protocols::jd2::Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin());

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	for ( ; it != job->output_string_string_pairs_end(); ++it ) {

		RowDataBaseOP key_data( new RowData<std::string>("data_key",it->first) );
		RowDataBaseOP value_data( new RowData<std::string>("data_value",it->second) );

		string_string_insert.add_row(utility::tools::make_vector(struct_id_data,key_data,value_data));
	}

	string_string_insert.write_to_database(db_session);

}

void
JobDataFeatures::load_string_string_data(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::pose::Pose &
){
	if ( !table_exists(db_session, "job_string_string_data") ) return;
	std::string statement_string =
		"SELECT\n"
		"\tdata_key,\n"
		"\tdata_value\n"
		"FROM\n"
		"\tjob_string_string_data\n"
		"WHERE\n"
		"\tjob_string_string_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while ( res.next() ) {
		std::string data_key, data_value;
		res >> data_key >> data_value;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string_string_pair(data_key, data_value);
	}
}

void JobDataFeatures::insert_string_real_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{

	InsertGenerator string_real_insert("job_string_real_data");
	string_real_insert.add_column("struct_id");
	string_real_insert.add_column("data_key");
	string_real_insert.add_column("data_value");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	protocols::jd2::Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin());

	for ( ; it != job->output_string_real_pairs_end(); ++it ) {

		RowDataBaseOP key_data( new RowData<std::string>("data_key",it->first) );
		RowDataBaseOP value_data( new RowData<core::Real>("data_value",it->second) );

		string_real_insert.add_row(utility::tools::make_vector(struct_id_data,key_data,value_data));
	}

	string_real_insert.write_to_database(db_session);
}

void
JobDataFeatures::load_string_real_data(
	utility::sql_database::sessionOP db_session,
	StructureID struct_id,
	core::pose::Pose &
){
	if ( !table_exists(db_session, "job_string_real_data") ) return;
	std::string statement_string =
		"SELECT\n"
		"\tdata_key,\n"
		"\tdata_value\n"
		"FROM\n"
		"\tjob_string_real_data\n"
		"WHERE\n"
		"\tjob_string_real_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while ( res.next() ) {
		std::string data_key;
		core::Real data_value;
		res >> data_key >> data_value;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string_real_pair(data_key, data_value);
	}
}


}
}
