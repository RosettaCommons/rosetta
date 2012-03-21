// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScoreTypeFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ScoreTypeFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>

// Platform Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <string>


namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::pose::Pose;
using core::graph::Graph;
using core::scoring::Energies;
using core::scoring::EnergyGraph;
using core::scoring::EnergyMap;
using core::scoring::EnergyEdge;
using core::scoring::getScoreFunction;
using core::scoring::hbond_sr_bb;
using core::scoring::hbond_lr_bb;
using core::scoring::n_score_types;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreType;
using core::scoring::ScoreFunctionOP;
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

ScoreTypeFeatures::ScoreTypeFeatures() :
	scfxn_(getScoreFunction())
{}

ScoreTypeFeatures::ScoreTypeFeatures(
	ScoreFunctionOP scfxn
) :
	scfxn_(scfxn)
{
	if ( scfxn_ == 0 ) {
		utility_exit_with_message( "ScoreTypeFeatures may not be constructed with a null-pointer ScoreFunctionOP" );
	}
}

ScoreTypeFeatures::ScoreTypeFeatures(
	ScoreTypeFeatures const & src
) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ScoreTypeFeatures::~ScoreTypeFeatures() {}

string
ScoreTypeFeatures::type_name() const { return "ScoreTypeFeatures"; }

string
ScoreTypeFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS score_types (\n"
			"	protocol_id INTEGER,\n"
			"	score_type_id INTEGER,\n"
			"	score_type_name TEXT,\n"
			"	PRIMARY KEY (protocol_id, score_type_id)\n"
			"	FOREIGN KEY (protocol_id)\n"
			"		REFERENCES protocols (protocol_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);\n"
			"\n";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS score_types (\n"
			"	protocol_id INTEGER REFERENCES protocols (protocol_id),\n"
			"	score_type_id INTEGER,\n"
			"	score_type_name TEXT,\n"
			"   PRIMARY KEY (protocol_id, score_type_id));\n"
			"\n";
	}else
	{
		return "";
	}

}

utility::vector1<std::string>
ScoreTypeFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}


Size
ScoreTypeFeatures::report_features(
	Size protocol_id,
	sessionOP db_session
){
	insert_score_type_rows(protocol_id, db_session);
	return 0;
}

void ScoreTypeFeatures::delete_record(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session){}

void
ScoreTypeFeatures::insert_score_type_rows(
	Size protocol_id,
	sessionOP db_session
) {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	std::string statement_string;
	if(db_mode == "sqlite3")
	{
		statement_string = "INSERT OR IGNORE INTO score_types VALUES (?,?,?);";
	}else if(db_mode == "mysql")
	{
		statement_string = "INSERT IGNORE INTO score_types VALUES (?,?,?);";
	}else
	{
		utility_exit_with_message("the database mode needs to be 'mysql' or 'sqlite3'");
	}

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));

		string const score_type( ScoreTypeManager::name_from_score_type(type) );
		stmt.bind(1,protocol_id);
		stmt.bind(2,score_type_id);
		stmt.bind(3,score_type);
		basic::database::safely_write_to_database(stmt);


	}
	//transact_guard.commit();
}


} // namesapce
} // namespace
