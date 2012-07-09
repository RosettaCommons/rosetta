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

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

//external headers
#include <cppdb/frontend.h>

#include <string>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

using basic::database::table_exists;

JobDataFeatures::JobDataFeatures() {}

JobDataFeatures::JobDataFeatures(JobDataFeatures const & ) {}

JobDataFeatures::~JobDataFeatures() {}

std::string JobDataFeatures::type_name() const
{
	return "JobDataFeatures";
}

void
JobDataFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;
	Column struct_id("struct_id",DbUUID(), false /*not null*/, false /*don't autoincrement*/);
	Column data_key("data_key",DbText(255));

	utility::vector1<Column> primary_columns;
	primary_columns.push_back(struct_id);
	primary_columns.push_back(data_key);
	PrimaryKey primary_key(primary_columns);

	Schema job_string_data("job_string_data",primary_key);
	job_string_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	job_string_data.write(db_session);

	Schema job_string_string_data("job_string_string_data",primary_key);
	job_string_string_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	job_string_string_data.add_column(Column("data_value",DbText()));

	job_string_string_data.write(db_session);

	Schema job_string_real_data("job_string_real_data",primary_key);
	job_string_real_data.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	job_string_real_data.add_column(Column("data_value",DbReal()));

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
		boost::uuids::uuid struct_id,
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
		boost::uuids::uuid struct_id,
		core::pose::Pose & pose
){
	load_string_data(db_session, struct_id, pose);
	load_string_string_data(db_session, struct_id, pose);
	load_string_real_data(db_session, struct_id, pose);
}

void JobDataFeatures::delete_record(
	boost::uuids::uuid struct_id,
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

void JobDataFeatures::insert_string_rows(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::Strings::const_iterator it(job->output_strings_begin());
	std::string statement_string = "INSERT INTO job_string_data (struct_id, data_key) VALUES (?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(; it != job->output_strings_end(); ++it)
	{
		stmt.bind(1,struct_id);
		stmt.bind(2,*it);
		basic::database::safely_write_to_database(stmt);
	}
}

void
JobDataFeatures::load_string_data(
		utility::sql_database::sessionOP db_session,
		boost::uuids::uuid struct_id,
		core::pose::Pose &
){
	if(!table_exists(db_session, "job_string_data")) return;

	std::string statement_string =
		"SELECT\n"
		"	data_key\n"
		"FROM\n"
		"	job_string_data\n"
		"WHERE\n"
		"	job_string_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next()){
		std::string data_key;
		res >> data_key;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string(data_key);
	}
}

void JobDataFeatures::insert_string_string_rows(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin());
	std::string statement_string = "INSERT INTO job_string_string_data (struct_id, data_key, data_value) VALUES (?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(; it != job->output_string_string_pairs_end();++it)
	{
		stmt.bind(1,struct_id);
		stmt.bind(2,it->first);
		stmt.bind(3,it->second);
		basic::database::safely_write_to_database(stmt);
	}
}

void
JobDataFeatures::load_string_string_data(
		utility::sql_database::sessionOP db_session,
		boost::uuids::uuid struct_id,
		core::pose::Pose &
){
	if(!table_exists(db_session, "job_string_string_data")) return;
	std::string statement_string =
		"SELECT\n"
		"	data_key,\n"
		"	data_value\n"
		"FROM\n"
		"	job_string_string_data\n"
		"WHERE\n"
		"	job_string_string_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next()){
		std::string data_key, data_value;
		res >> data_key >> data_value;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string_string_pair(data_key, data_value);
	}
}

void JobDataFeatures::insert_string_real_rows(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin());
	std::string statement_string = "INSERT INTO job_string_real_data (struct_id, data_key, data_value) VALUES (?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));


	for(; it != job->output_string_real_pairs_end();++it)
	{
		stmt.bind(1,struct_id);
		stmt.bind(2,it->first);
		stmt.bind(3,it->second);
		basic::database::safely_write_to_database(stmt);

	}
}

void
JobDataFeatures::load_string_real_data(
		utility::sql_database::sessionOP db_session,
		boost::uuids::uuid struct_id,
		core::pose::Pose &
){
	if(!table_exists(db_session, "job_string_real_data")) return;
	std::string statement_string =
		"SELECT\n"
		"	data_key,\n"
		"	data_value\n"
		"FROM\n"
		"	job_string_real_data\n"
		"WHERE\n"
		"	job_string_real_data.struct_id = ?;";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	while(res.next()){
		std::string data_key;
		core::Real data_value;
		res >> data_key >> data_value;
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string_real_pair(data_key, data_value);
	}
}


}
}
