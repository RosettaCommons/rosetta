// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProtocolFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/StructureFeatures.hh>

// Project Headers
#include <basic/database/sql_utils.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/functional/hash.hpp>

// C++
#include <string>
#include <sstream>
#include <iostream>




namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using basic::database::safely_prepare_statement;
using boost::hash_value;
using core::Size;
using core::io::silent::BinaryProteinSilentStruct;
using core::pose::Pose;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StructureFeatures::StructureFeatures(){}

StructureFeatures::StructureFeatures( StructureFeatures const & src) :
	FeaturesReporter()
{}

StructureFeatures::~StructureFeatures(){}

string
StructureFeatures::type_name() const { return "StructureFeatures"; }

string
StructureFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3")
		{
		return
			"CREATE TABLE IF NOT EXISTS structures (\n"
			"	struct_id INTEGER PRIMARY KEY AUTOINCREMENT,\n"
			"	protocol_id INTEGER,\n"
			"	tag TEXT,\n"
			"	input_tag TEXT,\n"
			"	UNIQUE (protocol_id, tag),"
			"	FOREIGN KEY (protocol_id)\n"
			"		REFERENCES protocols (protocol_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS structures (\n"
			"	struct_id BIGINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,\n"
			"	protocol_id INTEGER REFERENCES protocols(protocol_id),\n"
			"	tag VARCHAR(255),\n"
			"	input_tag VARCHAR(255),\n"
			"	UNIQUE (protocol_id, tag));";
			//"	FOREIGN KEY (protocol_id) REFERENCES protocols (protocol_id));";
	}else
	{
		return "";
	}
}

utility::vector1<std::string>
StructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}


//@details missing struct_id and input/output tags
Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size protocol_id,
	sessionOP db_session
){
	string const output_tag(protocols::jd2::JobDistributor::get_instance()->current_output_name());
	string const input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
	Size struct_id(
		report_features(pose, relevant_residues, protocol_id,
			db_session, output_tag, input_tag));
	return struct_id;
}

//@details missing struct_id and input/output tags
Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size protocol_id,
	sessionOP db_session,
	string const & tag,
	string const & input_tag
){
	string statement_string = "INSERT INTO structures VALUES (?,?,?,?);";
	statement stmt(safely_prepare_statement(statement_string,db_session));

	BinaryProteinSilentStruct silent_struct(pose, "");
	stringstream pose_string;
	silent_struct.print_conformation(pose_string);
	Size const struct_id = hash_value(pose_string.str());

	stmt.bind(1, struct_id);
	stmt.bind(2, protocol_id);
	stmt.bind(3, tag);
	stmt.bind(4, input_tag);
	basic::database::safely_write_to_database(stmt);

	return struct_id;
}


void StructureFeatures::delete_record(
	core::Size struct_id,
	utility::sql_database::sessionOP db_session
){

	std::string statement_string = "DELETE FROM structures WHERE struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(stmt);

}

void
StructureFeatures::load_into_pose(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	load_tag(db_session, struct_id, pose);
}

void
StructureFeatures::load_tag(
	sessionOP db_session,
	Size struct_id,
	Pose & pose) {

	std::string statement_string =
		"SELECT\n"
		"	tag\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.struct_id=?";

	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session))	;
	stmt.bind(1,struct_id);


	result res(basic::database::safely_read_from_database(stmt));
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with struct_id '"
			<< struct_id << "'." << endl;
		utility_exit_with_message(error_message.str());
	}
	string tag;
	res >> tag;
	tag_into_pose(pose,tag);
}

Size
StructureFeatures::get_struct_id(
	sessionOP db_session,
	string const & tag,
	core::Size const & protocol_id
){

	std::string statement_string =
		"SELECT\n"
		"	struct_id\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.tag=? AND structures.protocol_id=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,tag);
	stmt.bind(2,protocol_id);

	result res(basic::database::safely_read_from_database(stmt));
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with tag '"<<tag<<"'."<<endl;
		utility_exit_with_message(error_message.str());
	}
	long long struct_id;
	res >> struct_id;
	return static_cast<Size>(struct_id);
}

} // namesapce
} // namespace
