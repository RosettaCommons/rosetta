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
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// External Headers
#include <cppdb/frontend.h>

// C++
#include <string>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::pose::Pose;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StructureFeatures::StructureFeatures(){}

StructureFeatures::StructureFeatures( StructureFeatures const & ) :
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
			"	tag TEXT UNIQUE,\n"
			"	FOREIGN KEY (protocol_id)\n"
			"		REFERENCES protocols (protocol_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS structures (\n"
			"	struct_id INTEGER PRIMARY KEY AUTO_INCREMENT,\n"
			"	protocol_id INTEGER REFERENCES protocols(protocol_id),\n"
			"	tag VARCHAR(255) UNIQUE,\n"
			"	FOREIGN KEY (protocol_id) REFERENCES protocols (protocol_id));";
	}else
	{
		return "";
	}
}

Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size protocol_id,
	sessionOP db_session
){
	return report_features(
		pose, relevant_residues, protocol_id, db_session, find_tag(pose));
}

Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	Size protocol_id,
	sessionOP db_session
){
	return report_features(
		pose, relevant_residues, struct_id, protocol_id, db_session, find_tag(pose));
}


Size
StructureFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	Size protocol_id,
	sessionOP db_session,
	string const & tag
){
	statement stmt = (*db_session)
		<< "INSERT INTO structures VALUES (NULL,?,?);"
		<< protocol_id
		<< tag;
	stmt.exec();
	return stmt.last_insert_id();
}

Size
StructureFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	Size struct_id,
	Size protocol_id,
	sessionOP db_session,
	string const & tag
){
	statement stmt = (*db_session)
		<< "INSERT INTO structures VALUES (?,?,?);"
		<< struct_id
		<< protocol_id
		<< tag;
	stmt.exec();
	return stmt.last_insert_id();
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

	result res = (*db_session) <<
		"SELECT\n"
		"	tag\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.struct_id=?" << struct_id;
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
	string const & tag
){
	result res = (*db_session) <<
		"SELECT\n"
		"	struct_id\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.tag=?;" << tag;
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with tag '"<<tag<<"'."<<endl;
		utility_exit_with_message(error_message.str());
	}
	Size struct_id;
	res >> struct_id;
	return struct_id;
}

} // namesapce
} // namespace
