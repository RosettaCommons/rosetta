// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/features/JobDataFeatures.cc
/// @author Sam DeLuca

//unit headers
#include <protocols/features/JobDataFeatures.hh>


//platform headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/database/sql_utils.hh>

//external headers
#include <cppdb/frontend.h>

#include <string>

namespace protocols {
namespace features {

using basic::database::table_exists;

JobDataFeatures::JobDataFeatures() {}

JobDataFeatures::JobDataFeatures(JobDataFeatures const & src) {}

JobDataFeatures::~JobDataFeatures() {}

std::string JobDataFeatures::type_name() const
{
	return "JobDataFeatures";
}

std::string JobDataFeatures::schema() const
{
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
				"CREATE TABLE IF NOT EXISTS job_string_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key TEXT,\n"
				"	FOREIGN KEY (struct_id)\n"
				"	REFERENCES structures(struct_id)\n"
				"	DEFERRABLE INITIALLY DEFERRED,\n"
				"	PRIMARY KEY (struct_id,data_key));\n"
				"\n"
				"CREATE TABLE IF NOT EXISTS job_string_string_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key TEXT,\n"
				"	data_value TEXT,\n"
				"	FOREIGN KEY (struct_id)\n"
				"	REFERENCES structures(struct_id)\n"
				"	DEFERRABLE INITIALLY DEFERRED,\n"
				"	PRIMARY KEY (struct_id,data_key));\n"
				"\n"
				"CREATE TABLE IF NOT EXISTS job_string_real_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key TEXT,\n"
				"	data_value REAL,\n"
				"	FOREIGN KEY (struct_id)\n"
				"	REFERENCES structures(struct_id)\n"
				"	DEFERRABLE INITIALLY DEFERRED,\n"
				"	PRIMARY KEY (struct_id, data_key));";
	}else if(db_mode == "mysql")
	{
		return
				"CREATE TABLE IF NOT EXISTS job_string_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key VARCHAR(255),\n"
				"	FOREIGN KEY (struct_id)	REFERENCES structures(struct_id),\n"
				"	PRIMARY KEY (struct_id,data_key));\n"
				"\n"
				"CREATE TABLE IF NOT EXISTS job_string_string_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key VARCHAR(255),\n"
				"	data_value TEXT,\n"
				"	FOREIGN KEY (struct_id)	REFERENCES structures(struct_id),\n"
				"	PRIMARY KEY (struct_id,data_key));\n"
				"\n"
				"CREATE TABLE IF NOT EXISTS job_string_real_data (\n"
				"	struct_id INTEGER,\n"
				"	data_key VARCHAR(255),\n"
				"	data_value REAL,\n"
				"	FOREIGN KEY (struct_id)	REFERENCES structures(struct_id),\n"
				"	PRIMARY KEY (struct_id, data_key));";
	}else
	{
		return "";
	}

}

core::Size
JobDataFeatures::report_features(
		core::pose::Pose const & /*pose */,
		utility::vector1<bool> const & /*relevant_residues*/,
		core::Size struct_id,
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
		Size struct_id,
		core::pose::Pose & pose
){
	load_string_data(db_session, struct_id, pose);
	load_string_string_data(db_session, struct_id, pose);
	load_string_real_data(db_session, struct_id, pose);
}

void JobDataFeatures::delete_record(
	core::Size struct_id,
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

void JobDataFeatures::insert_string_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::Strings::const_iterator it(job->output_strings_begin());
	std::string statement_string = "INSERT INTO job_string_data VALUES (?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	for(; it != job->output_strings_end(); ++it)
	{
		stmt.bind(2,*it);
		basic::database::safely_write_to_database(stmt);
	}
}

void
JobDataFeatures::load_string_data(
		utility::sql_database::sessionOP db_session,
		Size struct_id,
		core::pose::Pose & pose
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

void JobDataFeatures::insert_string_string_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin());
	std::string statement_string = "INSERT INTO job_string_string_data VALUES (?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	for(; it != job->output_string_string_pairs_end();++it)
	{
		stmt.bind(2,it->first);
		stmt.bind(3,it->second);
		basic::database::safely_write_to_database(stmt);
	}
}

void
JobDataFeatures::load_string_string_data(
		utility::sql_database::sessionOP db_session,
		Size struct_id,
		core::pose::Pose & pose
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

void JobDataFeatures::insert_string_real_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const
{
	protocols::jd2::Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin());
	std::string statement_string = "INSERT INTO job_string_real_data VALUES (?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);

	for(; it != job->output_string_real_pairs_end();++it)
	{
		stmt.bind(2,it->first);
		stmt.bind(3,it->second);
		basic::database::safely_write_to_database(stmt);

	}
}

void
JobDataFeatures::load_string_real_data(
		utility::sql_database::sessionOP db_session,
		Size struct_id,
		core::pose::Pose & pose
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

