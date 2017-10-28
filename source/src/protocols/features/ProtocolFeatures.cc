// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProtocolFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProtocolFeatures.hh>

// Project Headers
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <core/types.hh>
#include <utility/version.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

//Basic Headers
#include <basic/Tracer.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <sstream>

static THREAD_LOCAL basic::Tracer TR( "protocols.features.ProtocolFeatures" );

namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using basic::options::OptionKeys::parser::protocol;
using basic::options::option;
using core::Size;
using core::pose::Pose;
using utility::io::izstream;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowData;
using basic::database::insert_statement_generator::RowDataBaseOP;

ProtocolFeatures::ProtocolFeatures(){}

ProtocolFeatures::ProtocolFeatures( ProtocolFeatures const & ) : utility::pointer::ReferenceCount() {}

ProtocolFeatures::~ProtocolFeatures()= default;

string
ProtocolFeatures::type_name() const { return "ProtocolFeatures"; }

void
ProtocolFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session, core::Size protocol_id) const{
	using namespace basic::database::schema_generator;

	//if protocol id is set, don't autoincrement
	if ( protocol_id ) {

		Column protocol_id("protocol_id", DbDataTypeOP( new DbInteger() ));
		Schema protocols("protocols", PrimaryKey(protocol_id));

		protocols.add_column( Column("specified_options", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("command_line", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("svn_url", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("svn_version", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("script", DbDataTypeOP( new DbText() )) );
		protocols.write(db_session);
	} else {
		//if protocol id is not set, don't autoincrement

		Column protocol_id("protocol_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, true /*autoincrement*/);
		Schema protocols("protocols", PrimaryKey(protocol_id));

		protocols.add_column( Column("specified_options", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("command_line", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("svn_url", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("svn_version", DbDataTypeOP( new DbText() )) );
		protocols.add_column( Column("script", DbDataTypeOP( new DbText() )) );
		protocols.write(db_session);
	}
}

utility::vector1<std::string>
ProtocolFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	return dependencies;
}


string
ProtocolFeatures::indices() const {
	return "";
}

Size
ProtocolFeatures::report_features(
	Size protocol_id,
	sessionOP db_session
){

	string const & command_line( option.get_argv() );

	stringstream option_stream;
	option_stream << basic::options::option;
	string const & specified_options( option_stream.str() );

	string const & svn_url( utility::Version::url() );
	string const & svn_version( utility::Version::version() );

	bool using_rosetta_scripts( basic::options::option[ protocol ].active() );
	string script = "";
	if ( using_rosetta_scripts ) {
		string const & script_fname( basic::options::option[ protocol ] );
		stringstream script_buf;
		script_buf << utility::io::izstream( script_fname.c_str() ).rdbuf();
		script = script_buf.str();
	}

	cppdb::statement insert_statement;
	if ( protocol_id ) {
		//Check to make sure a protocol with the same id doesn't already exist (This should only happen if you are manually setting the protocol_id
		//through the options system or the ReportToDB tag in RosettaScripts
		std::string statement_string =
			"SELECT\n"
			"\tcount(*)\n"
			"FROM\n"
			"\tprotocols\n"
			"WHERE\n"
			"\tprotocols.protocol_id = ?;";
		cppdb::statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,protocol_id);

		TR << "Checking for existing protocol entry with given id" << std::endl;
		cppdb::result res(basic::database::safely_read_from_database(stmt));
		if ( res.next() ) {
			core::Size selected = 0;
			res >> selected;
			if ( selected != 0 ) {
				return protocol_id;
			}
		}

		TR << "Writing to protocols table with given protocol id: " << protocol_id << std::endl;
		std::string insert_string("INSERT INTO protocols (protocol_id, specified_options, command_line, svn_url, svn_version, script) VALUES (?,?,?,?,?,?);");
		insert_statement = basic::database::safely_prepare_statement(insert_string,db_session);
		insert_statement.bind(1,protocol_id);
		insert_statement.bind(2,command_line);
		insert_statement.bind(3,specified_options);
		insert_statement.bind(4,svn_url);
		insert_statement.bind(5,svn_version);
		insert_statement.bind(6,script);

	} else {
		TR << "No protocol ID, generating one automagically" << std::endl;
		std::string insert_string("INSERT INTO protocols (command_line, specified_options, svn_url, svn_version, script) VALUES (?,?,?,?,?);");
		insert_statement = basic::database::safely_prepare_statement(insert_string,db_session);
		insert_statement.bind(1,command_line);
		insert_statement.bind(2,specified_options);
		insert_statement.bind(3,svn_url);
		insert_statement.bind(4,svn_version);
		insert_statement.bind(5,script);
	}

	//Protocol features doesn't use safely_write_to_database due to special handling of thrown cppdb exceptions
	//basic::database::safely_write_to_database(insert_statement);
	insert_statement.exec();
	if ( protocol_id ) {
		return protocol_id;
	} else {
		core::Size new_protocol_id = insert_statement.sequence_last("protocols_protocol_id_seq");
		return new_protocol_id;
	}
}

} // features namesapce
} // protocols namespace
