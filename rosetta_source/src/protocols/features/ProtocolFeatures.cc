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
#include <protocols/features/ProtocolFeatures.hh>

// Project Headers
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/svn_version.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using basic::options::OptionKeys::parser::protocol;
using basic::options::option;
using core::Size;
using core::minirosetta_svn_url;
using core::minirosetta_svn_version;
using core::pose::Pose;
using utility::io::izstream;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

ProtocolFeatures::ProtocolFeatures(){}

ProtocolFeatures::ProtocolFeatures( ProtocolFeatures const & ) :
	FeaturesReporter()
{}

ProtocolFeatures::~ProtocolFeatures(){}

string
ProtocolFeatures::type_name() const { return "ProtocolFeatures"; }

string
ProtocolFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	bool protocol_id_mode = basic::options::option[basic::options::OptionKeys::out::database_protocol_id].user();
	if(db_mode == "sqlite3")
	{
		if(protocol_id_mode)
		{
			return
				"CREATE TABLE IF NOT EXISTS protocols (\n"
				"	protocol_id INTEGER PRIMARY KEY UNIQUE,\n"
				"	command_line TEXT,\n"
				"	specified_options TEXT,\n"
				"	svn_url TEXT,\n"
				"	svn_version TEXT,\n"
				"	script TEXT);";
		}else
		{
			//default behavior
			return
				"CREATE TABLE IF NOT EXISTS protocols (\n"
				"	protocol_id INTEGER PRIMARY KEY AUTOINCREMENT,\n"
				"	command_line TEXT,\n"
				"	specified_options TEXT,\n"
				"	svn_url TEXT,\n"
				"	svn_version TEXT,\n"
				"	script TEXT);";
		}

	}else if(db_mode == "mysql")
	{
		if(protocol_id_mode)
		{
			return
				"CREATE TABLE IF NOT EXISTS protocols (\n"
				"	protocol_id INTEGER PRIMARY KEY UNIQUE,\n"
				"	command_line TEXT,\n"
				"	specified_options TEXT,\n"
				"	svn_url TEXT,\n"
				"	svn_version TEXT,\n"
				"	script TEXT);";
		}else
		{
			//default behavior
			return
				"CREATE TABLE IF NOT EXISTS protocols (\n"
				"	protocol_id INTEGER PRIMARY KEY AUTO_INCREMENT,\n"
				"	command_line TEXT,\n"
				"	specified_options TEXT,\n"
				"	svn_url TEXT,\n"
				"	svn_version TEXT,\n"
				"	script TEXT);";
		}

	}else
	{
		return "";
	}
}

string
ProtocolFeatures::indices() const {
	return "";
}

Size
ProtocolFeatures::report_features(
	Pose const & /*pose*/,
	vector1< bool > const & /*relevant_residues*/,
	Size protocol_id,
	sessionOP db_session
){

	string const command_line( option.get_argv() );

	stringstream option_stream;
	option_stream << basic::options::option;
	string const specified_options( option_stream.str() );

	string const svn_url( minirosetta_svn_url() );
	string const svn_version( minirosetta_svn_version() );

	bool using_rosetta_scripts( basic::options::option[ protocol ].active() );
	string script = "";
	if ( using_rosetta_scripts){
		string const script_fname( basic::options::option[ protocol ] );
		stringstream script_buf;
		script_buf << utility::io::izstream( script_fname.c_str() ).rdbuf();
		script = script_buf.str();
	}

	//if -out:database_protocol_id is specified we need to make sure the protocol hasn't already been specified

	cppdb::statement stmt = (*db_session) <<
		"SELECT\n"
		"	count(*)\n"
		"FROM\n"
		"	protocols\n"
		"WHERE\n"
		"	protocols.protocol_id == ?;" << protocol_id;

	cppdb::result res(basic::database::safely_read_from_database(stmt));
	if(res.next())
	{
		core::Size selected = 0;
		res >> selected;
		if(selected != 0)
		{
			return protocol_id;
		}
	}


	if(protocol_id){
		stmt = (*db_session)
			<< "INSERT INTO protocols VALUES (?,?,?,?,?,?);"
			<< protocol_id;
	} else {
		stmt = (*db_session)
			<< "INSERT INTO protocols VALUES (NULL,?,?,?,?,?);";
	}
	stmt
		<< command_line
		<< specified_options
		<< svn_url
		<< svn_version
		<< script;
	basic::database::safely_write_to_database(stmt);
	return stmt.sequence_last("");

}

} // features namesapce
} // protocols namespace
